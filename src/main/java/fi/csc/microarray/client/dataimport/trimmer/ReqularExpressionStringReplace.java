package fi.csc.microarray.client.dataimport.trimmer;

import java.util.regex.Pattern;
import java.util.regex.PatternSyntaxException;

public class ReqularExpressionStringReplace extends DataTrimmingOperation{
	
	private String regexp;
	private String newString;

	public ReqularExpressionStringReplace(String regexp, String newString, int columnIndex) throws PatternSyntaxException {
		super(columnIndex);
		
		// Test the regular expression
		Pattern.compile(regexp);
		
		this.regexp = regexp;
		this.newString = newString;
	}
	
	@Override
	public String doTrimming(String stringToTrim) {
		return stringToTrim.replaceAll(regexp, newString);
	}
}
