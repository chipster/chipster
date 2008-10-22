package fi.csc.microarray.util;

import java.util.Iterator;

public class Deseparator implements Iterable<String[]>, Iterator<String[]> {
	
	private AdvancedStringTokenizer ast;
	private int batchSize;
	private boolean atEnd = false;
	private String separator;
	
	public Deseparator(String separator, AdvancedStringTokenizer ast, int batchSize) {
		this.separator = separator;
		this.ast = ast;
		this.batchSize = batchSize;
	}
	
	public Iterator<String[]> iterator() {
		return this;
	}
	
	public boolean hasNext() {
		return !atEnd;
	}
	
	public String[] next() {
		assert(!atEnd);
		
		// read in a batch
		String[] batch = new String[batchSize];
		for (int i = 0; i < batch.length; i++) {				
			batch[i] = ast.next();
		}
		
		// check for end conditions
		if (batch[batchSize-1].endsWith(separator)) {
			batch[batchSize-1] = batch[batchSize-1].substring(0, batch[batchSize-1].length()-1);
			
		} else if (ast.hasNext() && ast.peek().equals(separator)) {
			ast.next();
			
		} else {
			atEnd = true;
		}
		return batch;
	}
	
	public void remove() {
		throw new UnsupportedOperationException();			
	}		
}
