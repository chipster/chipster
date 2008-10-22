package fi.csc.microarray.util;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.StringTokenizer;

public class AdvancedStringTokenizer implements Iterable<String>, Iterator<String> {

	private List<String> tokens = new LinkedList<String>();
	int index = -1;
	private int inputStringLength;

	public AdvancedStringTokenizer(String string) {
		this(string, true, false, " \t\n\r\f");		
	}
	
	public AdvancedStringTokenizer(String string, boolean trim, boolean groupHyphenated) {
		this(string, trim, groupHyphenated, " \t\n\r\f");				
	}

	/**
	 * 
	 * @param string string to tokenise
	 * @param trim remove leading and trailing space from tokens
	 * @param groupHyphenated group text inside hyphens into one token
	 * @param delim characters that delimit tokens
	 */
	public AdvancedStringTokenizer(String string, boolean trim, boolean groupHyphenated, String delim) {
		StringTokenizer st = new StringTokenizer(string, delim);
		this.inputStringLength = string.length();
		
		while (st.hasMoreTokens()) {
			String token = trim ? st.nextToken().trim() : st.nextToken();
			tokens.add(token);
		}
		
		if (groupHyphenated) {
			List<String> originalTokens = tokens; 
			tokens = new LinkedList<String>();
			for (int i = 0; i < originalTokens.size(); i++) {
				if (!originalTokens.get(i).startsWith("\"")) {
					tokens.add(originalTokens.get(i));
				} else {
					String groupedToken = "";
					while (!originalTokens.get(i).endsWith("\"")) {						
						groupedToken += (originalTokens.get(i) + " ");
						i++;
					}
					groupedToken += originalTokens.get(i);
					tokens.add(groupedToken.substring(1, groupedToken.length() - 1));
				}
			}
		}
	}

	public Iterator<String> iterator() {
		return this; // just to make this Iterator also Iterable and usable in for:each
	}

	public boolean hasNext() {
		return (index+1) < tokens.size();
	}

	public String next() {
		if (!hasNext()) {
			throw new RuntimeException("tokeniser error: no tokens left (there were " + tokens.size() + " tokens in input string of length " + inputStringLength + ")");
		}
		return tokens.get(++index);
	}

	public String peek() {
		return tokens.get(index + 1);
	}

	public String current() {
		return tokens.get(index);
	}

	/**
	 * Not supported.
	 * 
	 * @throws UnsupportedOperationException always
	 */
	public void remove() {
		throw new UnsupportedOperationException();
	}
}
