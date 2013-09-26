package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * This class represents a single line of Vcf file.
 * 
 * @author klemela
 */
public class TsvLine extends FileLine implements SelectionText  {

	private Region region;
	private String[] headers;
	private String[] values;

	public Region getRegion() {
		return region;
	}
	public void setRegion(Region region) {
		this.region = region;
	}
	public String[] getHeaders() {
		return headers;
	}
	public void setHeaders(String[] headers) {
		this.headers = headers;
	}
	public String[] getValues() {
		return values;
	}
	public void setValues(String[] values) {
		this.values = values;
	}

	@Override
	public String getText() {				
	
		String text = "TSV feature info\n";
		
		for (int i = 0; i < headers.length && i < values.length; i++) {
			text += headers[i] + "\t" + values[i] + "\n";
		}
		
		return text;				
	}
}