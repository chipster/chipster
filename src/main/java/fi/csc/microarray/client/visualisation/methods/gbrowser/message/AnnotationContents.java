package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.Collection;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

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
		public String file;

		public Row(String species, String version, String content, String file) {
			this.species = species;
			this.version = version;		
			this.file = file;
			
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

	public void parseFrom(InputStream contentsStream) throws IOException {

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
			rows.add(new Row(splitted[0], splitted[1], splitted[2], splitted[3]));
		}

	}


	public void add(Row row) {
		rows.add(row);
	}
	

	public void write() throws IOException {
		contentsFile.delete();
		Writer writer = null;
		try {
			writer = new FileWriter(contentsFile, true);
			writer.write(FILE_ID + "\n");
			for (Row row : rows) {
				writer.write(row.species + "\t" + row.version + "\t" + row.content + "\t" + row.file + "\n");
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

	public Collection<Genome> getGenomes() {
		Set<Genome> genomes = new LinkedHashSet<Genome>();
		for (Row row : rows) {
			genomes.add(row.getGenome());
		}
		return genomes;
	}
}
