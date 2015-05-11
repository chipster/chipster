package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.net.URLConnection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;
import org.yaml.snakeyaml.Yaml;
import org.yaml.snakeyaml.constructor.Constructor;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.KeyAndTrustManager;

/**
 * 
 * Class for parsing annotation information from 
 * the folder structure of the public file broker urls.
 * 
 * @author Petri Klemelä, Taavi Hupponen
 * 
 */
public class AnnotationManager {
	
	private static final String ANNOTATIONS = "genomebrowser";

	private static final Logger logger = Logger.getLogger(AnnotationManager.class);
	
	//Location parts of the external genome browser urls are replaced with these strings in the genomes.yaml file
	public static final String CHR_LOCATION = "[CHR]";
	public static final String START_LOCATION = "[START]";
	public static final String END_LOCATION = "[END]";

	private File localAnnotationsRoot;

	private LinkedList<GenomeAnnotation> annotations = new LinkedList<GenomeAnnotation>();
	private GBrowser browser;

	private List<URL> remoteAnnotationFiles;

	private HashMap<String, String> displaySpecies = new HashMap<String, String>();
	private HashMap<String, String> displayVersions = new HashMap<String, String>();

	private List<Genome> genomes = new LinkedList<>();

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
		 * @return
		 */
		public DataUrl getUrl() {
			if (checkLocalFile(this)) {
				String fileName = getLocalFileName(this);
				File localFile = new File(localAnnotationsRoot, fileName);

				try {
					return new DataUrl(localFile);
				} catch (MalformedURLException e) {
					logger.warn("generating url for local file " + localFile + "failed");
				}
			}
			return new DataUrl(this.url, url.getPath());
		}

		public long getContentLength() {
			return this.contentLength;
		}

		public Genome getGenome() {
			return new Genome(species, version);
		}
		
		public String toString() {
			return getClass().getSimpleName() + " " + getGenome() + " " + getUrl();
		}
	}

	public static class Genome implements Comparable<Genome> {
		public String speciesId;
		public String versionId;
		public String displaySpecies;
		public String displayVersion;
		public String sortId;

		public Genome(String species, String version) {
			this.speciesId = species;
			this.versionId = version;
		}

		@Override
		public String toString() {
						
			return getSpecies() + " " + getVersion();
		}
		
		private String getSpecies() {
			if (displaySpecies != null) {
				return displaySpecies;
			} else {
				return speciesId;
			}
		}
		
		private String getVersion() {
			if (displayVersion != null) {
				return displayVersion;
			} else {
				return versionId;
			}
		}

		@Override
		public boolean equals(Object o) {
			if (o instanceof Genome) {
				Genome other = (Genome) o;
				return compareTo(other) == 0;
			}
			return false;
		}
		
		@Override
		public int hashCode() {
			//Eclipse generated implementation
			final int prime = 31;
			int result = 1;
			result = prime * result
					+ ((speciesId == null) ? 0 : speciesId.hashCode());
			result = prime * result
					+ ((versionId == null) ? 0 : versionId.hashCode());
			return result;
		}
		

		@Override
		public int compareTo(Genome other) {
			int speciesComparison = speciesId.compareTo(other.speciesId);
			int versionComparison;

			Integer ensemblVersion =  parseEnsemblVersion(versionId);
			Integer otherEnsemblVersion =  parseEnsemblVersion(other.versionId);
			
			if (ensemblVersion != null && otherEnsemblVersion != null) {
				versionComparison = -ensemblVersion.compareTo(otherEnsemblVersion);
			} else {
				// alphabetical order if parsing fails
				versionComparison = -versionId.compareTo(other.versionId);
			}
			
			if (speciesComparison != 0){
				return speciesComparison;
			} else {
				return versionComparison;
			}
		}

		private Integer parseEnsemblVersion(String versionId) {			
			int period = versionId.lastIndexOf(".");
			
			if (period >= 0) {
				String ensemblString = versionId.substring(period + 1);
				try {
					return Integer.parseInt(ensemblString);
				} catch (NumberFormatException e) {
					// just return null
				}
			}			
			return null;
		}
	}
	
	public enum AnnotationType {
		CYTOBANDS("Cytoband"), 
		GTF_TABIX("Transcript"), GTF_TABIX_INDEX("Transcript index"), REPEAT("Repeat"), REPEAT_INDEX("Repeat index"),
		REFERENCE("Reference sequence", false), REFERENCE_INDEX("Reference sequence index"), SNP("ENSEMBL SNP"), SEARCH_INDEX("Search index"), 
		ENSEMBL_BROWSER_URL("Ensembl", false), UCSC_BROWSER_URL("UCSC", false), GENOME_INFO("Annotation definition", true), GTF("Gene annotation", false);

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
		this.remoteAnnotationFiles = browser.getRemoteAnnotationFiles();
		this.localAnnotationsRoot = browser.getLocalAnnotationDir();
		
		if (this.remoteAnnotationFiles != null) {
			interpretAnnotationFiles(remoteAnnotationFiles);
		} else {
			// remote annotations weren't available, try local files
			logger.info("trying to use local annotations");
			
			//parse directories
			LinkedList<URL> urls = new LinkedList<URL>();
			
			for (File file : getLocalFiles(localAnnotationsRoot)) {
				urls.add(file.toURI().toURL());
			}
			
			interpretAnnotationFiles(urls);
		}
		
		if (genomes.isEmpty()) {
			throw new GBrowserException("Cannot access genome browser annotations from the server or from the local cache.");
		}

		removeUnnecessaryFiles(localAnnotationsRoot);
	}

	private boolean interpretAnnotationFiles(List<URL> files) throws URISyntaxException, MalformedURLException {
		
		boolean isSuccess = false;
		
		for (URL file: files) {
			
			String path = file.toURI().getPath();//decode url, e.g. convert %20 to space
			
			path = path.substring(path.indexOf(ANNOTATIONS) + ANNOTATIONS.length());			
			
			String[] directories = path.split("/");
			
			//this is a proper directory-defined annotation only if there are exactly two sub-directories (species and version)
			//(split creates two additional Strings due to the preceding and trailing slashes)
			if (directories.length == 4) {
				
				isSuccess = true;
				
				//String empty = directories[0];
				String species = directories[1];
				String version = directories[2];
				String fileName = directories[3];
				
				AnnotationType annotationType = null;
				
				if (fileName.endsWith(".search.tsv")) {
					annotationType = AnnotationType.SEARCH_INDEX;
					
				} else if (fileName.endsWith(".tabix.gtf.gz")) {
					annotationType = AnnotationType.GTF_TABIX;
					
				} else if (fileName.endsWith(".tabix.gtf.gz.tbi")) {
					annotationType = AnnotationType.GTF_TABIX_INDEX;
					
				} else if (fileName.endsWith(".gtf")) {
					annotationType = AnnotationType.GTF;
					
				} else if (fileName.endsWith("cytoband-chr.txt")) {
					annotationType = AnnotationType.CYTOBANDS;
					
				} else if (fileName.endsWith("repeat-tabix.bed.gz")) {
					annotationType = AnnotationType.REPEAT;
					
				} else if (fileName.endsWith("repeat-tabix.bed.gz.tbi")) {
					annotationType = AnnotationType.REPEAT_INDEX;
					
				} else if (fileName.endsWith(".fa")) {
					annotationType = AnnotationType.REFERENCE;
					
				} else if (fileName.endsWith(".fa.fai")) {
					annotationType = AnnotationType.REFERENCE_INDEX;
					
				} else if (fileName.endsWith(".yaml")) {

					annotationType = AnnotationType.GENOME_INFO;
					parseGenomeInfo(file, species, version);										
				}  									
				
				if (annotationType != null) {
					
					GenomeAnnotation annotation = new GenomeAnnotation(species, version, annotationType.getId(), null, file, -1);
					addAnnotation(annotation);
					
					//If there was no sortId in GenomeInfo, then generate one.
					//These are sorted alphabetically (after defined sortStrings and before contents2 genomes) 
					Genome genome = new Genome(species, version);
					if (!genomes.contains(genome)) {
						genome.sortId =  "\uFFFE";
						genomes.add(genome);
					}
				}
			}			
		}
		
		return isSuccess;
	}

	private void parseGenomeInfo(URL file, String species, String version) {
			
		try {
			Yaml yaml = new Yaml(new Constructor(GenomeInfo.class));
			URLConnection connection = file.openConnection();
			KeyAndTrustManager.configureForChipsterCertificate(connection);
			GenomeInfo info = (GenomeInfo) yaml.load(connection.getInputStream());
			
			if (info != null) {
				
				if (info.getSpecies() != null) {
					displaySpecies.put(species, info.getSpecies());
				}
				
				if (info.getVersion() != null) {
					displayVersions.put(version, info.getVersion());
				}
				
				if (info.getEnsembl() != null) {

					URL ensemblUrl = info.getEnsembl();
					GenomeAnnotation ensembl = new GenomeAnnotation(
							species, version, AnnotationType.ENSEMBL_BROWSER_URL.getId(), null, ensemblUrl, -1);
					addAnnotation(ensembl);
				}

				if (info.getBrowserUrl() != null) {

					URL ucscUrl = info.getBrowserUrl();
					GenomeAnnotation ucsc = new GenomeAnnotation(
							species, version, AnnotationType.UCSC_BROWSER_URL.getId(), null, ucscUrl, -1);			
					addAnnotation(ucsc);
				}
				
				if (info.getSortId() != null) {
					Genome genome = new Genome(species, version);
					genome.sortId = info.getSortId();

					//Remove possible genome with generated sortid
					genomes.remove(genome);
					genomes.add(genome);
				}
			}

		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void removeUnnecessaryFiles(File localAnnotationsRoot) throws IOException {

		List<File> allFiles = getLocalFiles(localAnnotationsRoot);

		for (File file : allFiles) {
		
			if (!this.contains(file)) {
				removeFile(file);
			}
		}
		
		removeEmptyDirectories(localAnnotationsRoot);
	}
	
	private void removeEmptyDirectories(File path) {
		
		for (File file : path.listFiles()) {
			
			if (file.isDirectory()) {
				removeEmptyDirectories(file);
				
				if(file.listFiles().length == 0) {
					try {
						this.removeFile(file);
					} catch (IOException e) {
						logger.error(e.getMessage(), e);
					}
				}
			}
		}
	}

	private void removeFile(File fileToRemove) throws IOException {
		//Just one more check, in case something is horribly wrong 
		if (fileToRemove.getCanonicalPath().contains(".chipster")) {
			fileToRemove.delete();
		} else {
			logger.error("Attempt to remove file '" + fileToRemove + "' was prevented because " +
					"the file is not inside folder '.chipster'");
		}
	}

	private boolean contains(File file) throws IOException {
				
		for (GenomeAnnotation annotation : annotations) {
			if (annotation.url != null) {
				
				String fileName = getLocalFileName(annotation);

				if (file.getCanonicalPath().endsWith(fileName)) {
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

		for (GenomeAnnotation annotation : annotations) {
						
			int index = genomes.indexOf(annotation.getGenome());
			
			if(index < 0) {
				logger.error("All annotation genomes are not in the genome list: " + annotation.getGenome());
			} else {
				Genome genome = genomes.get(index);
				
				if (displaySpecies.containsKey(genome.speciesId)) {
					genome.displaySpecies = displaySpecies.get(genome.speciesId);
				}
				
				if (displayVersions.containsKey(genome.versionId)) {
					genome.displayVersion = displayVersions.get(genome.versionId);
				}				
			}								
		}
		
		Collections.sort(genomes, new Comparator<Genome>() {
			@Override
			public int compare(Genome o1, Genome o2) {
				int sortIdComparison = o1.sortId.compareTo(o2.sortId);
				if (sortIdComparison != 0) {
					return sortIdComparison;
				} else {
					return o1.compareTo(o2);
				}				
			}			
		});
		
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
								downloadAnnotationFile(annotation);
								
							} catch (IOException e) {
								throw new RuntimeException(e);
							}
						}
					}
				}
			}
		});
	}

	private void downloadAnnotationFile(GenomeAnnotation annotation) throws IOException {
		
		String fileName = getLocalFileName(annotation);
				
		File localFile = new File(this.localAnnotationsRoot, fileName);
		
		localFile.getParentFile().mkdirs();
		
		InputStream in = null;
		try {
			URLConnection connection = annotation.url.openConnection();
			KeyAndTrustManager.configureForChipsterCertificate(connection);
			in = connection.getInputStream();
			IOUtils.copy(in, localFile);
		} finally {
			IOUtils.closeIfPossible(in);
		}
	}

	private String getLocalFileName(GenomeAnnotation annotation) {
		
		Genome genome = annotation.getGenome();		
		String species = genome.speciesId;
		String version = genome.versionId;
		
//		//replace everything that is not a word character (a-z in any case, 0-9 or _)
//		species = species.replaceAll("[^\\w]"," ").trim();
//		version = version.replaceAll("[^\\w]"," ").trim();
//		
//		//remove repeating spaces
//		species = species.replaceAll("\\s+","-");
//		version = version.replaceAll("\\s+","-");
	
		String fileName = IOUtils.getFilenameWithoutPath(annotation.url);
		
		return File.separator +  species + File.separator + version + File.separator + fileName;
	}

	private boolean checkLocalFile(GenomeAnnotation annotation) {
		if (annotation.url != null) {
			String fileName = getLocalFileName(annotation);
			File localFile = new File(this.localAnnotationsRoot, fileName);
			
			if (localFile.exists() ) {
				
				//check that the content lengths match, if we know it (length is not negative)
				if (annotation.getContentLength() < 0 || localFile.length() == annotation.getContentLength()) {
					return true;
				} else {
					//There was a file with same name than the annotation, but the size differs.
					//Propably it's just some old version of the annotation, so let's just remove it to 
					//make space for downloading a new version.
					try {
						removeFile(localFile);
					} catch (IOException e) {
						logger.error(e.getMessage(), e);
					}
				}
			}
		}
		return false;
	}

	public void addAnnotation(GenomeAnnotation annotation) {
		annotations.add(annotation);
	}
	
	public List<File> getLocalFiles(File path) throws IOException {
		
		List<File> fileList = new LinkedList<File>();
		
		addFilesRecursively(fileList, path);
				
		return fileList;
	}
	
	private void addFilesRecursively(List<File> files, File path) throws IOException {
		
		for (File file : path.listFiles()) {
			
			if (file.isDirectory()) {
				addFilesRecursively(files, file);
				
			} else {

				files.add(file);
			}
		}		
	}
}
