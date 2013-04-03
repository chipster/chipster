package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;


public abstract class TsvLineParser implements LineParser {
	
	protected String[] values;
	
	public Long getLong(int column) {
		String string = values[column];
		return new Long(string);		
	}
	
	public Float getFloat(int column) {
		String string = values[column];
		return new Float(string);		
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
	
	public boolean isContentLine() {
		return values != null;
	}

	public abstract String getHeaderStart();
}
