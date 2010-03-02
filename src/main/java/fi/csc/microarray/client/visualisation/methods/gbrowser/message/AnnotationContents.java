package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

/**
 * 
 * Class for saving and loading contents file. 
 * 
 * @author klemela
 *
 */
public class AnnotationContents {

	private final File contentsFile = new File("contents.txt");

	private final String FILE_ID = "CHIPSTER ANNOTATION CONTENTS FILE VERSION 1";



	public void add(Row row) {
		appendToFile(row.species + "\t" + row.version + "\t" + row.content + "\t" + row.file);
	}

	public List<Row> getContents(InputStream contentsStream) {

		List<Row> list = new ArrayList<Row>();

		try {

			BufferedReader reader = new BufferedReader(new InputStreamReader(contentsStream));

			reader.readLine().equals(FILE_ID);

			String line;
			while((line = reader.readLine()) != null) {

				String[] splitted = line.split("\t");

				list.add(new Row(splitted[0], splitted[1], splitted[2], splitted[3]));
			}

		} catch (FileNotFoundException e) {

			e.printStackTrace();
		} catch (IOException e) {

			e.printStackTrace();
		}

		return list;
	}

	private void appendToFile(String str) {

		try {

			Writer writer = new FileWriter(contentsFile, true);
			writer.write(str + "\n");

			writer.flush();
			writer.close();

		} catch (IOException e) {

			e.printStackTrace();
		} 
	}

	public void clear() {
		contentsFile.delete();
		appendToFile(FILE_ID);
	}
	
	/**
	 * Data structure for Contents class
	 * 
	 * @author klemela
	 *
	 */
	public class Row {
		public Row(String species, String version, String content,
				String file) {
			this.species = species;
			this.version = version;
			this.content = content;
			this.file = file;
		}
		public String species;
		public String version;
		public String content;
		public String file;
	}

}

