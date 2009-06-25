package fi.csc.microarray.analyser.bsh;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.LookaheadLineReader;

/**
 * Utilities to be used by Java and BeanShell jobs.
 * 
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 *
 */
// TODO partly duplicates functionality in Feature API, but do we want to start using DataManager and DataBeans here?
public class JavaJobUtils {

	/**
	 * Get the gene names from GENELIST type of input file.
	 * 
	 * If there is only one column, use it. If there are many columns,
	 * search for "symbol" column, if not found, search for " " and "identifier" columns.
	 * 
	 * @param file
	 * @return the gene names, empty String[] if the file is empty or gene names are not found
	 * @throws IOException
	 */
	public static String[] getGeneNames(File file) throws IOException {
		return getColumns(file, new String[] {"symbol", " ", "identifier"});
	}

	public static String[] getProbes(File file) throws IOException {
		return getColumns(file, new String[] {" "});
	}

	public static String[] getColumns(File file, String[] lookForColumns) throws IOException {

		List<String> names = new LinkedList<String>();
		LookaheadLineReader reader = new LookaheadLineReader(new BufferedReader(new FileReader(file)));
		
		try {
			String line;
			String[] columns;
			line = reader.readLine();

			// empty file
			if (line == null) {
				return new String[0];
			}
			
			// split first line
			columns = line.split("\t");
			
			// check that column count is not different from header due to invisible row name column
			if (reader.peekLine().split("\t").length == (columns.length + 1)) {
				String[] newColumns = new String[columns.length + 1];
				System.arraycopy(columns, 0, newColumns, 1, columns.length);
				newColumns[0] = " "; // must be space, empty names are not allowed
				columns = newColumns;
			}
		
			int nameColumnIndex = -1;
			for (String lookForColumn : lookForColumns) {
				nameColumnIndex = Arrays.asList(columns).indexOf(lookForColumn);
				if (nameColumnIndex != -1) {
					break; // found, skip the rest
				}
			}
			
			if (nameColumnIndex == -1) {
				throw new IOException("data does not contain any of the required columns: " + Arrays.asList(lookForColumns));
			}
			
			// read in the rest of the file
			for (line = reader.readLine(); line != null; line = reader.readLine()) {
				names.add(line.split("\t")[nameColumnIndex]);
			}
			
		} finally {
			IOUtils.closeIfPossible(reader.getReader());
		}

		return names.toArray(new String[names.size()]);
	}
}
