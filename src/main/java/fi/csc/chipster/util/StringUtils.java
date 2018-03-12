package fi.csc.chipster.util;

public class StringUtils {

	public static String capitalizeFirstLetter(String original) {
	    if (original == null || original.length() == 0) {
	        return original;
	    }
	    return original.substring(0, 1).toUpperCase() + original.substring(1);
	}
	
	public static String removeLastLine(String s) {
		int i = s.lastIndexOf("\n");
		if (i != -1) {
			return s.substring(0, i);
		} else {
			return s;
		}
	}
	
}
