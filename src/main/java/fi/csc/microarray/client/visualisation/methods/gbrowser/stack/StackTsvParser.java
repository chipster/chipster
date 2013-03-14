package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;


public abstract class StackTsvParser implements Parser {
	
	protected String[] values;
	
	public Long getLong(int column) {
		String string = values[column];
		return new Long(string);		
	}
	
	public String getString(int column) {
		return values[column];
	}
	
	@Override
	public boolean setLine(String line) {
		if (line.startsWith(getHeaderStart())) {
			this.values = null;
			return false;
		} else {
			this.values = line.split("\t"); 
			return true; 
		}
	}
	
	public boolean isContentLinle() {
		return values != null;
	}

	public abstract String getHeaderStart();
}
