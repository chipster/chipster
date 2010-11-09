package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.Writer;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.util.IOUtils;

/**
 * 
 * Class for saving and loading contents file. The contents file describes what annotation data files we have
 * available at the annotation repository.
 * 
 * @author Petri Klemelä
 * 
 */
public class AnnotationContents {
	private static final String CONTENTS_FILE = "contents.txt";
	private static final String ANNOTATIONS_PATH = "annotations";

	private static final Logger logger = Logger.getLogger(AnnotationContents.class);

	
	private URL remoteAnnotationsRoot;
	private File localAnnotationsRoot;
	
	private final File contentsFile = new File("contents.txt");

	private final String FILE_ID = "CHIPSTER ANNOTATION CONTENTS FILE VERSION 1";

	private LinkedList<Row> rows = new LinkedList<Row>();

	/**
	 * Data structure for Contents class
	 * 
	 */
	public class Row {
		
		public String species;
		public String version;
		public Content content;
		public URL url;

		public Row(String species, String version, String content, URL url) {
			this.species = species;
			this.version = version;		
			this.url = url;
			
			for (Content type : Content.values()) {
				if (type.id.equals(content)) {
					this.content = type;
				}
			}
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
		
		@Override 
		public int hashCode() {
			return species.hashCode();
		}
	}
	
	public enum Content { 
		CYTOBANDS ("Cytobands"), 
		TRANSCRIPTS ("ENSEMBL Transcripts"),
		GENES ("ENSEMBL Genes"),
		MIRNA ("ENSEMBL miRNA Genes"),
		REFERENCE ("Reference sequence"),
		SNP ("ENSEMBL SNP");
		
		String id;
		
		Content(String id) {
			this.id = id;
		}

		public String getId() {
			return id;
		}
	}

	
	
	public void initialize() throws Exception {

		// get annotation locations
		this.remoteAnnotationsRoot = getRemoteAnnotationsUrl();
		this.localAnnotationsRoot = DirectoryLayout.getInstance().getLocalAnnotationDir();
		
		// try to parse the remote contents file
		URL remoteContents = IOUtils.createURL(remoteAnnotationsRoot, CONTENTS_FILE);
		boolean remoteContentsOk;
		InputStream remoteContentsStream = null;
		try {
			remoteContentsStream = remoteContents.openStream();
			parseFrom(remoteContentsStream);
			remoteContentsOk = true;
		} catch (Exception e) {
			remoteContentsOk = false;
		} finally {
			IOUtils.closeIfPossible(remoteContentsStream);
		}
		
		// if everything went well, also make a local copy of contents file
		// it will be used when working offline
		File localContents = new File(localAnnotationsRoot, CONTENTS_FILE);
		if (remoteContentsOk) {
			logger.info("using remote annotation contents file");
			OutputStream localContentsStream = null;;
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
				throw e;
			} finally {
				IOUtils.closeIfPossible(localContentsStream);
			}
		}
	}

	
	
	public void write() throws IOException {
		contentsFile.delete();
		Writer writer = null;
		try {
			writer = new FileWriter(contentsFile, true);
			writer.write(FILE_ID + "\n");
			for (Row row : rows) {
				writer.write(row.species + "\t" + row.version + "\t" + row.content + "\t" + row.url.getFile() + "\n");
			}
		} finally {
			IOUtils.closeIfPossible(writer);
		}
	}

	public List<Row> getRows() {
		return rows;
	}
	
	public Row getRow(Genome genome, Content content) {

		for (Row row : rows) {
			if (row.getGenome().equals(genome) && row.content == content) {
				return row;
			}
		}
		return null;
	}

	public List<Genome> getGenomes() {
		List<Genome> genomes = new LinkedList<Genome>();
		for (Row row : rows) {
			if (!genomes.contains(row.getGenome())) {
				genomes.add(row.getGenome());
			}
		}
		return genomes;
	}
	
	/**
	 * Returns local if there are no annotations.
	 * 
	 * 
	 * TODO What if there are no annotations except for reference?
	 * @param genome
	 * @return
	 */
	public boolean hasLocalAnnotations(Genome genome) {
		for (Content c : Content.values()) {
			if (!c.equals(Content.REFERENCE)) {
				Row annotation = getRow(genome, c);
				if (annotation != null && !IOUtils.isLocalFileURL(annotation.url)) {
					return false;
				}
			}
		}
		return true;
	}

	/**
	 * Returns local if there is no reference.
	 * 
	 * @param genome
	 * @return
	 */
	public boolean hasLocalReference(Genome genome) {
		Row reference = getRow(genome, Content.REFERENCE);
		if (reference != null && !IOUtils.isLocalFileURL(reference.url)) {
			return false;
		} else {
			return true;
		}
	}

	public void downloadAnnotations(Genome genome) throws IOException {
		for (Content c : Content.values()) {
			if (!c.equals(Content.REFERENCE)) {
				Row annotation = getRow(genome, c);
				if (annotation != null && !checkLocalFile(annotation)) {
					downloadAnnotationFile(annotation.url);
				}
			}
		}
	}

	private void downloadAnnotationFile(URL sourceUrl) throws IOException {
		String fileName = sourceUrl.getFile();
		File localFile = new File(this.localAnnotationsRoot, fileName);
		InputStream in = null;
		try  {
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
	private boolean checkLocalFile(Row annotation) {
		if (annotation != null && IOUtils.isLocalFileURL(annotation.url)) {
			return true;
		} else {
			return false;
		}
	}

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
			
			// use local file if it exists
			URL url;
			String fileName = splitted[3];
			File localFile = new File(localAnnotationsRoot, fileName);
			// FIXME add more checks, size and checksum maybe
			if (localFile.exists()) {
				url = localFile.toURI().toURL();
			} else {
				url = IOUtils.createURL(remoteAnnotationsRoot, fileName);
			}
			rows.add(new Row(splitted[0], splitted[1], splitted[2], url));
		}
	
	}



	private URL getRemoteAnnotationsUrl() throws Exception {
		FileBrokerClient fileBroker = Session.getSession().getServiceAccessor()
				.getFileBrokerClient();
		URL annotationsUrl = new URL(fileBroker.getPublicUrl() + "/"
				+ ANNOTATIONS_PATH);
		return annotationsUrl;
	}

	public static void main(String[] args) throws IOException {
		URL url = new URL("http://chipster-filebroker.csc.fi:8090/public/annotations/Homo_sapiens.NCBI36.54_transcripts.tsv");
		System.out.println(url.openConnection().getContentLength());
	}
	
}
