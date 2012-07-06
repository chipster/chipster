package fi.csc.microarray.manager.web.util;

public class StringUtils {
	
	public static String getHumanReadable(long longValue) {
		
		String stringValue;
		
    	if (longValue >= 1000000000000l) {
    		stringValue = "" + (longValue / 1000000000000l) + " TB";
    		
    	} else if (longValue >= 1000000000) {
    		stringValue = "" + (longValue / 1000000000) + " GB";
    		
    	} else if (longValue >= 1000000) {
    		stringValue = "" + (longValue / 1000000) + " MB";
    		
    	} else {
    		stringValue = "" + (longValue / 1000) + " kB";
    	}
    	
    	return stringValue;
	}
}
