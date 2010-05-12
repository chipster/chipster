package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.util.IOUtils;

/**
 * 
 * Class for saving and loading contents file.
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
		public String content;
		public String file;

		public Row(String species, String version, String content, String file) {
			this.species = species;
			this.version = version;
			this.content = content;
			this.file = file;
		}
	}

	public void parseFrom(InputStream contentsStream) {

		List<Row> list = new ArrayList<Row>();

		try {

			BufferedReader reader = new BufferedReader(new InputStreamReader(contentsStream));

			reader.readLine().equals(FILE_ID);

			String line;
			while ((line = reader.readLine()) != null) {
				String[] splitted = line.split("\t");
				list.add(new Row(splitted[0], splitted[1], splitted[2], splitted[3]));
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace(); // TODO fix exception handling
			
		} catch (IOException e) {
			e.printStackTrace(); // TODO fix exception handling
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


	public LinkedHashSet<String> getGenomes() {
		LinkedHashSet<String> genomes = new LinkedHashSet<String>();
		for (Row row : rows) {
			genomes.add(row.species + " " + row.version);
		}
		return genomes;
	}
}
