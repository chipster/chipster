package fi.csc.microarray.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.util.LinkedList;

public class LookaheadLineReader {
	
	private BufferedReader reader;
	private LinkedList<String> buffer = new LinkedList<String>(); 
	
	public LookaheadLineReader(BufferedReader reader) {
		this.reader = reader;	
	}
	
	public BufferedReader getReader() {
		return reader;
	}
	
	public String readLine() throws IOException {
		if (buffer.isEmpty()) {
			return reader.readLine();
		} else {
			return buffer.poll();
		}
	}
	
	public String peekLine() throws IOException {
		return peekLine(1);
	}
	
	public String peekLine(int lookahead) throws IOException {
		// fill buffer if needed
		while (buffer.size() < lookahead) {
			buffer.offer(reader.readLine());
		}
		return buffer.get(lookahead-1);
		
	}

	public void read(int chars) throws IOException {
		String line = peekLine();
		line = line.substring(chars);
		if (line.length() == 0) {
			readLine(); // remove first line because it would become empty
		} else {
			// replace line with shortened version
			buffer.remove(0);
			buffer.add(0, line); 
		}
	}

}
