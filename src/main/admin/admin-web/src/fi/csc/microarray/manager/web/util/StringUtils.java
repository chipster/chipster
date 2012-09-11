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
	
	/**
	 * Converts minutes to hours and minutes separated by colon, i.e. H:MM
	 * 
	 * @param minutes
	 * @return hours and minutes separated by colon, i.e. H:MM
	 */
	public static String formatMinutes(int minutes) {
		int hours = minutes / 60;
		int minuteRemainder = minutes % 60;
		
		return String.format("%d:%02d", hours, minuteRemainder);
	}
}
