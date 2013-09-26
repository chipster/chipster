package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;


public abstract class AbstractTsvLineParser implements LineParser {
	
	protected String[] values;

	public Integer getInteger(int column) {
		String string = values[column];
		return new Integer(string);		
	}
	
	public Long getLong(int column) {
		String string = values[column];
		return new Long(string);		
	}
	
	public Float getFloat(int column) {
		String string = values[column];
		try {
			return new Float(string);
		} catch (NumberFormatException e) {
			if ("inf".equals(string.toLowerCase())) {
				return Float.POSITIVE_INFINITY;
			} else if ("-inf".equals(string.toLowerCase())) {
				return Float.NEGATIVE_INFINITY;
			}
			return Float.NaN;
		}
	}
	
	public String getString(int column) {
		return values[column];
	}
	
	@Override
	public boolean setLine(String line) {
		if (getHeaderStart() != null && line.startsWith(getHeaderStart())) {
			this.values = null;
			return false;
		} else {
			this.values = line.split("\t"); 
			return true; 
		}
	}
	
	@Override
	public boolean isContentLine() {
		return values != null;
	}

	public String getHeaderStart() {
		return null;
	}
}
