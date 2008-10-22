package fi.csc.microarray.client.dataimport.trimmer;

public class ConditionalStringReplace extends DataTrimmingOperation{

	private double lowerLimit;
	private double upperLimit;
	private String replacement;
	private boolean acceptEqualsFirst;
	private boolean acceptEqualsLast;

	/**
	 * Creates conditional string replace operation. The acceptEquals parameters 
	 * might be a bit weird at the first sight. They works like this: <br><br>
	 * <ul>
	 * <li><code>acceptEqualsFirst: <b>false</b>, acceptEqualsLast: <b>false</b>: </code> lowerLimit < x < upperLimit</li>
	 * <li><code>acceptEqualsFirst: <b>true</b>, acceptEqualsLast: <b>false</b>: </code> lowerLimit <= x < upperLimit</li>
	 * <li><code>acceptEqualsFirst: <b>false</b>, acceptEqualsLast: <b>true</b>: </code> lowerLimit < x <= upperLimit</li>
	 * <li><code>acceptEqualsFirst: <b>true</b>, acceptEqualsLast: <b>true</b>: </code> lowerLimit <= x <= upperLimit</li>
	 * 
	 * @param lowerLimit
	 * @param upperLimit
	 * @param acceptEqualsFirst < or =< for first condition
	 * @param acceptEqualsLast < or =< for last condition
	 * @param replacement
	 * @param columnIndex
	 */
	public ConditionalStringReplace(double lowerLimit, double upperLimit, boolean acceptEqualsFirst, boolean acceptEqualsLast, String replacement, int columnIndex) {
		super(columnIndex);
		this.lowerLimit = lowerLimit;
		this.upperLimit = upperLimit;
		this.acceptEqualsFirst = acceptEqualsFirst;
		this.acceptEqualsLast = acceptEqualsLast;
		this.replacement = replacement;
	}

	@Override
	public String doTrimming(String stringToTrim) {
		double value;
		try{
			value = Double.parseDouble(stringToTrim.trim());
		} catch (NumberFormatException nfe) {
			return stringToTrim;
		}
		
		if(acceptEqualsFirst == acceptEqualsLast){
			if(acceptEqualsFirst){
				// Both true
				if(lowerLimit <= value && value <= upperLimit){
					return replacement;
				} else {
					return stringToTrim;
				}
			} else {
				// Both false
				if(lowerLimit < value && value < upperLimit){
					return replacement;
				} else {
					return stringToTrim;
				}
			}
		} else {
			if(acceptEqualsFirst){
				// First true
				if(lowerLimit <= value && value < upperLimit){
					return replacement;
				} else {
					return stringToTrim;
				}
			} else {
				// Last true
				if(lowerLimit < value && value <= upperLimit){
					return replacement;
				} else {
					return stringToTrim;
				}
			}
		}
	}
	
}
