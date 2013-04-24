package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.util.IOUtils;

/**
 * 
 * Class for saving and loading annotation contents file. The contents file describes what
 * annotation data files we have available at the annotation repository.
 * 
 * @author Petri Klemelä, Taavi Hupponen
 * 
 */
public class AnnotationManager {
	private static final String CONTENTS_FILE = "contents2.txt";

	private static final Logger logger = Logger.getLogger(AnnotationManager.class);
	
	//Location parts of the external genome browser urls are replaced with these strings in the contents file
	public static final String CHR_LOCATION = "[CHR]";
	public static final String START_LOCATION = "[START]";
	public static final String END_LOCATION = "[END]";

	private URL remoteAnnotationsRoot;
	private File localAnnotationsRoot;

	private final String FILE_ID = "CHIPSTER ANNOTATION CONTENTS FILE VERSION 2";
	private final String CHR_UNSPECIFIED =  "*";

	private LinkedList<GenomeAnnotation> annotations = new LinkedList<GenomeAnnotation>();
	private GBrowser browser;

	/**
	 * Model for single annotation file.
	 * 
	 */
	public class GenomeAnnotation {

		public String species;
		public String version;
		public AnnotationType type;
		public Chromosome chr;

		private URL url;
		private long contentLength;

		public GenomeAnnotation(String species, String version, String annotationType, Chromosome chr, URL url, long contentLength) {
			this.species = species;
			this.version = version;
			this.chr = chr;
			this.contentLength = contentLength;
			this.url = url;

			for (AnnotationType aType : AnnotationType.values()) {
				if (aType.id.equals(annotationType)) {
					this.type = aType;
				}
			}
		}

		/**
		 * Return url pointing to a local file if the local file ok.
		 * 
		 * TODO what to do if offline and contents.txt content length and local
		 * file content length mismatch
		 * 
		 * @return
		 */
		public URL getUrl() {
			if (checkLocalFile(this)) {
				String fileName = IOUtils.getFilenameWithoutPath(this.url);
				File localFile = new File(localAnnotationsRoot, fileName);
				URL newUrl;
				try {
					newUrl = localFile.toURI().toURL();
				} catch (MalformedURLException e) {
					logger.warn("generating url for local file " + localFile + "failed");
					return this.url;
				}
				return newUrl;
			}
			return this.url;
		}

		public long getContentLength() {
			return this.contentLength;
		}

		public Genome getGenome() {
			return new Genome(species, version);
		}
	}

	public class Genome {
		public String species;
		public String version;

		public Genome(String species, String version) {
			this.species = species;
			this.version = version;
		}

		@Override
		public String toString() {
			return species + " " + version;
		}

		@Override
		public boolean equals(Object o) {
			if (o instanceof Genome) {
				Genome other = (Genome) o;
				return species.equals(other.species) && version.equals(other.version);
			}
			return false;
		}

		// FIXME not good
		@Override
		public int hashCode() {
			return species.hashCode();
		}
	}
	
	public enum AnnotationType {
		CYTOBANDS("Cytoband"), 
		GTF_TABIX("Transcript"), GTF_TABIX_INDEX("Transcript index"), REPEAT("Repeat"), REPEAT_INDEX("Repeat index"),
		REFERENCE("Reference sequence", false), REFERENCE_INDEX("Reference sequence index"), SNP("ENSEMBL SNP"), GENE_CHRS("Gene name"), 
		ENSEMBL_BROWSER_URL("Ensembl", false), UCSC_BROWSER_URL("UCSC", false);

		private String id;
		private boolean clientCacheable;

		AnnotationType(String id) {
			this(id, true);
		}
		
		AnnotationType(String id, boolean clientCacheable) {
			this.id = id;
			this.clientCacheable = clientCacheable;
		}

		public String getId() {
			return id;
		}

		public boolean isClientCacheable() {
			return clientCacheable;
		}
	}

	public AnnotationManager(GBrowser browser) {
		this.browser = browser;
	}

	/**
	 * Get and parse the contents.txt, which describes available annotations.
	 * 
	 * @throws Exception
	 */
	public void initialize() throws Exception {

		// get annotation locations
		this.remoteAnnotationsRoot = browser.getRemoteAnnotationsUrl();
		this.localAnnotationsRoot = browser.getLocalAnnotationDir();

		// try to parse the remote contents file
		boolean remoteContentsOk = false;
		InputStream remoteContentsStream = null;
		URL remoteContents = null;
		if (this.remoteAnnotationsRoot != null) {
			remoteContents = IOUtils.createURL(remoteAnnotationsRoot, CONTENTS_FILE);
			try {
				remoteContentsStream = remoteContents.openStream();
				parseFrom(remoteContentsStream);
				remoteContentsOk = true;
			} catch (Exception e) {
				remoteContentsOk = false;
			} finally {
				IOUtils.closeIfPossible(remoteContentsStream);
			}
		}

		// if everything went well, also make a local copy of contents file
		// it will be used when working offline
		File localContents = new File(localAnnotationsRoot, CONTENTS_FILE);
		if (remoteContentsOk) {
			logger.info("using remote annotation contents file");
			OutputStream localContentsStream = null;

			try {
				remoteContentsStream = remoteContents.openStream();
				localContentsStream = new FileOutputStream(localContents);
				IOUtils.copy(remoteContentsStream, localContentsStream);
			} catch (Exception e) {
				logger.warn("could not make a local copy of contents file", e);
			} finally {
				IOUtils.closeIfPossible(remoteContentsStream);
				IOUtils.closeIfPossible(localContentsStream);
			}
		}

		// remote contents could not be loaded, try local contents file
		else {
			logger.info("trying to use local annotation contents file");
			InputStream localContentsStream = null;
			try {
				localContentsStream = new BufferedInputStream(new FileInputStream(localContents));
				parseFrom(localContentsStream);
			} catch (Exception e) {
				// also local contents file failed
				throw new Exception("Cannot access genome browser annotations from the server or from the local cache.", e);
			} finally {
				IOUtils.closeIfPossible(localContentsStream);
			}
		}

		//if (remoteContentsOk) {
		removeUnnecessaryFiles(localAnnotationsRoot);
		//}
	}

	private void removeUnnecessaryFiles(File localAnnotationsRoot) {
		File annotationFolder = localAnnotationsRoot;

		String[] allFiles = annotationFolder.list();

		for (String file : allFiles) {
			if (!this.contains(file)) {
				File fileToRemove = new File(localAnnotationsRoot, file);
				removeFile(fileToRemove);
			}
		}
	}
	
	private void removeFile(File fileToRemove) {
		//Just one more check, in case something is horribly wrong 
		if (fileToRemove.getPath().contains(".chipster")) {
			fileToRemove.delete();
		} 
	}

	private boolean contains(String file) {
		
		if (file.equals(CONTENTS_FILE)) {
			return true;
		}
		
		for (GenomeAnnotation annotation : annotations) {
			if (annotation.url != null) {
				String path = annotation.url.getPath();
				String fileName = path.substring(path.lastIndexOf("/") + 1);

				if (fileName.equals(file)) {
					return true;
				}
			}
		}
		return false;
	}

	public List<GenomeAnnotation> getAnnotations() {
		return annotations;
	}

	public List<GenomeAnnotation> getAnnotations(Genome genome, AnnotationType annotationType) {
		List<GenomeAnnotation> filteredAnnotations = new LinkedList<GenomeAnnotation>();
		for (GenomeAnnotation annotation : annotations) {

			if (annotation.getGenome().equals(genome) && annotation.type == annotationType) {
				filteredAnnotations.add(annotation);
			}
		}
		return filteredAnnotations;
	}

	public GenomeAnnotation getAnnotation(Genome genome, AnnotationType annotationType) {

		List<GenomeAnnotation> filteredList =  getAnnotations(genome, annotationType);

		if (filteredList.size() > 0) {
			return filteredList.get(0);
		} else {
			return null;
		}
	}


	public List<Genome> getGenomes() {
		List<Genome> genomes = new LinkedList<Genome>();
		for (GenomeAnnotation annotation : annotations) {
			if (!genomes.contains(annotation.getGenome())) {
				genomes.add(annotation.getGenome());
			}
		}
		return genomes;
	}

	/**
	 * Returns local if there are no annotations. In such a case there is
	 * nothing to be downloaded.
	 * 
	 * @param genome
	 * @return
	 */
	public boolean hasLocalAnnotations(Genome genome) {
		for (AnnotationType c : AnnotationType.values()) {
			if (c.isClientCacheable()) {
				GenomeAnnotation annotation = getAnnotation(genome, c);
				if (annotation != null && !checkLocalFile(annotation)) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * Download annotations for a genome. Downloading happens as a blocking
	 * task.
	 * 
	 * @param genome
	 * @throws IOException
	 */
	public void downloadAnnotations(final Genome genome) throws IOException {
		browser.runBlockingTask("downloading annotations", new Runnable() {

			@Override
			public void run() {
				for (AnnotationType c : AnnotationType.values()) {
					if (c.isClientCacheable()) {
						GenomeAnnotation annotation = getAnnotation(genome, c);
						if (annotation != null && !checkLocalFile(annotation)) {

							// don't use getUrl() here because we need the
							// remote url
							try {
								downloadAnnotationFile(annotation.url);
							} catch (IOException e) {
								throw new RuntimeException(e);
							}
						}
					}
				}
			}
		});
	}

	private void downloadAnnotationFile(URL sourceUrl) throws IOException {
		String fileName = sourceUrl.getPath().substring(sourceUrl.getPath().lastIndexOf('/') + 1);
		File localFile = new File(this.localAnnotationsRoot, fileName);
		InputStream in = null;
		try {
			in = sourceUrl.openStream();
			IOUtils.copy(in, localFile);
		} finally {
			IOUtils.closeIfPossible(in);
		}
	}

	/**
	 * TODO add check for file size and or checksum
	 * 
	 */
	private boolean checkLocalFile(GenomeAnnotation annotation) {
		if (annotation.url != null) {
			String fileName = IOUtils.getFilenameWithoutPath(annotation.url);
			File localFile = new File(this.localAnnotationsRoot, fileName);
			if (localFile.exists() ) {
				if (localFile.length() == annotation.getContentLength()) {
					return true;
				} else {
					//There was a file with same name than the annotation, but the size differs.
					//Propably it's just some old version of the annotation, so let's just remove it to 
					//make space for downloading a new version.
					removeFile(localFile);
				}
			}
		}
		return false;
	}

	/**
	 * Parse contents file.
	 * 
	 * @param contentsStream
	 * @throws IOException
	 */
	private void parseFrom(InputStream contentsStream) throws IOException {

		BufferedReader reader = new BufferedReader(new InputStreamReader(contentsStream));

		if (!reader.readLine().equals(FILE_ID)) {
			throw new IllegalArgumentException("annotation stream does not start with " + FILE_ID);
		}

		String line;
		while ((line = reader.readLine()) != null) {
			if (line.trim().equals("")) {
				continue;
			}
			String[] splitted = line.split("\t");

			// Try to always store the remote url even if a local file exists.
			// Existence of the local is checked later every time it is needed.
			URL url;
			String fileName = splitted[4];
			
			if ("".equals(fileName)) {
				url = null;
			} else if (fileName.startsWith("http://")) {
				//Not a real filename, but a full url
				url = new URL(fileName);
			} else {
				url = IOUtils.createURL(remoteAnnotationsRoot != null ? remoteAnnotationsRoot : new URL("file://"), fileName);
			}
			
			long contentLength = Long.parseLong(splitted[5]);

			Chromosome chr = null;
			if (!splitted[3].equals(CHR_UNSPECIFIED)) {
				chr = new Chromosome(splitted[3]);
			}

			addAnnotation(new GenomeAnnotation(splitted[0], splitted[1], splitted[2], chr, url, contentLength));
		}
	}

	public void addAnnotation(GenomeAnnotation annotation) {
		annotations.add(annotation);
	}
}
