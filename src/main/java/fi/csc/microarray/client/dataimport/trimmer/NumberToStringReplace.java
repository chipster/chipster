package fi.csc.microarray.client.dataimport.trimmer;

/**
 * Replace operation to replace numerical flag value to string value
 * 
 * @author mkoski
 *
 */
public class NumberToStringReplace extends DataTrimmingOperation {

	private Double oldNumber;
	private String newString;

	public NumberToStringReplace(Double oldNumber, String newString, int columnIndex) {
		super(columnIndex);
		this.oldNumber = oldNumber;
		this.newString = newString;
	}
	
	@Override
	public String doTrimming(String stringToTrim) {
		try{
		Double numberToReplace = Double.valueOf(stringToTrim);
			if(numberToReplace.equals(oldNumber)){
				return newString;
			} else {
				return stringToTrim;
			}
		} catch(NumberFormatException nfe){
			return stringToTrim;
		}
	}
}
