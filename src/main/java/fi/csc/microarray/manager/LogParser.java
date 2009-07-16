package fi.csc.microarray.manager;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;


/**
 * TODO Need way to close inputStream??
 * 
 * 
 * @author hupponen
 *
 */
public class LogParser {


	public Iterable<HashMap<String, String>> parse(File file) throws IOException {
		LogIterable iterable = new LogIterable(new FileInputStream(file)); 
		
		return iterable;
	}
	

	private static HashMap<String, String> parseEntry(String entry) {
		HashMap<String, String> result = new HashMap<String, String>();
		
		String[] tokens = entry.split(";");
		
		// check for old log files with , as delimiter, remove check
		// after changing delimiter
		if (tokens.length == 7) {
			result.put("id", tokens[0]);
			result.put("operation", tokens[0]);
			result.put("operation", tokens[1]);
			result.put("status", tokens[2]);
			result.put("username", tokens[3]);
			result.put("starttime", tokens[4]);
			result.put("endtime", tokens[5]);
			result.put("compHost", tokens[6]);
		} else {
			throw new RuntimeException("Illegal number of tokens: " + tokens.length);
		}
		return result;
	}

	
	
	@SuppressWarnings("unused")
	private static HashMap<String, String> parseOldEntry(String entry) {
		HashMap<String, String> result = new HashMap<String, String>();
		
		String[] tokens = entry.split(",");
		
		// check for old log files with , as delimiter, remove check
		// after changing delimiter
		if (tokens.length == 5) {
			result.put("operation", tokens[0]);
			result.put("status", tokens[1]);
			result.put("username", tokens[2]);
			result.put("starttime", tokens[3]);
			result.put("endtime", tokens[4]);
		} else if (tokens.length == 6 && tokens[0].equals("Annotation/Agilent")) {
			result.put("operation", tokens[0] + "," + tokens[1]);
			result.put("status", tokens[2]);
			result.put("username", tokens[3]);
			result.put("starttime", tokens[4]);
			result.put("endtime", tokens[5]);
		}
		else {
			throw new RuntimeException("Illegal number of tokens: " + tokens.length);
		}
		return result;
	}

	
	
	



	private class LogIterable implements Iterable<HashMap<String, String>> {
		
		private class LineIterator implements Iterator<HashMap<String, String>> {

			private BufferedReader reader;
			private String nextLine;
			
			public LineIterator(InputStream input) throws IOException {
				this.reader = new BufferedReader(new InputStreamReader(input));
				nextLine = reader.readLine();
				if (nextLine == null) {
					reader.close();
				}
			}
			
			public boolean hasNext() {
				return nextLine != null;
			}

			public HashMap<String, String> next() {
				if (nextLine == null) {
					return null;
				}
				
				HashMap<String, String> entry = LogParser.parseEntry(nextLine);
				try {
					nextLine = reader.readLine();
					if (nextLine == null) {
						reader.close();
					}
				} catch (IOException e) {
					throw new RuntimeException(e);
				}
				return entry;
			}

			public void remove() {
				throw new UnsupportedOperationException();
			}
		}

		
		
		
		private LineIterator iterator;
		
		public LogIterable(InputStream inputStream) throws IOException { 
			this.iterator = new LineIterator(inputStream);
		}
		
		public Iterator<HashMap<String, String>> iterator() {
			return iterator;
		}
		
	}

	
	
	
	
	public static void main(String[] args) throws IOException {
		LogParser parser = new LogParser();
		Iterable<HashMap<String, String>> entries = parser.parse(new File("logfile.log"));
		for (HashMap<String, String> entry: entries) {
			System.out.println(entry.get("username"));
		}
	}
	
}
