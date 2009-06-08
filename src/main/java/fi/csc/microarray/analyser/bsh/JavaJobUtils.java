package fi.csc.microarray.analyser.bsh;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.util.IOUtils;

/**
 * Utilities to be used by Java and BeanShell jobs.
 * 
 * 
 * @author hupponen
 *
 */
public class JavaJobUtils {

	/**
	 * Get the gene names from GENELIST type of input file.
	 * 
	 * If there is only one column, use it. If there are many columns,
	 * search for "symbol" column, if not found, search for " " column.
	 * 
	 * @param file
	 * @return the gene names, empty String[] if the file is empty or gene names are not found
	 * @throws IOException
	 */
	public static String[] getGeneNames(File file) throws IOException {

		List<String> names = new LinkedList<String>();
		BufferedReader reader = new BufferedReader(new FileReader(file));
		try {
			String line;
			String[] columns;
			int nameColumnIndex = 0;
			line = reader.readLine();

			// empty file
			if (line == null) {
				return new String[0];
			}
			
			// split first line
			columns = line.split("\t");
			
			// only one column
			if (columns.length == 1) {
				nameColumnIndex = 0;
			} 
				
			// many columns, find the gene names column
			else if ( columns.length > 1) {
				// search for symbol column
				nameColumnIndex = Arrays.asList(columns).indexOf("symbol");

				// symbol not found, search for "" column
				if (nameColumnIndex < 0) {
					nameColumnIndex = Arrays.asList(columns).indexOf(" ");

					// "" not found, give up
					if (nameColumnIndex < 0) {
						return new String[0];
					}
				}
			}				
			
			
			// read in the rest of the file
			for (line = reader.readLine(); line != null; line = reader.readLine()) {
				names.add(line.split("\t")[nameColumnIndex]);
			}
		} finally {
			IOUtils.closeIfPossible(reader);
		}

		return names.toArray(new String[names.size()]);
	}

}
