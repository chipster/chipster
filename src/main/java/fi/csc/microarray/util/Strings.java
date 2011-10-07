package fi.csc.microarray.util;

import java.util.StringTokenizer;

public class Strings {
		
	/**
	 * Combines strings into one string by placing a delimeter between them.
	 * Works with objects by calling toString() for them.
	 * If there are no strings to combine, empty string is returned.
	 */
	public static String delimit(Iterable<? extends Object> objects, String delimeter) {
		String result = "";
		boolean first = true;
		for (Object object : objects) {
			if (first) {
				result = object.toString();
				first = false;
			} else {
				result += (delimeter + object);
			}
		}
		return result;		
	}
	
	public static String repeat(String string, int times) {
		String result = "";
		for ( int c = 0; c < times; c++ ) {
		      result += string;
		}
		return result;
	}

	public static String[] repeatToArray(String string, int times) {
		String[] result = new String[times];
		for ( int c = 0; c < times; c++ ) {
		      result[c] = string;
		}
		return result;
	}

	public static String indent(String string, int spaces) {
		StringTokenizer st = new StringTokenizer(string, "\n");
		String result = "";
		while (st.hasMoreTokens()) {
			result += (repeat(" ", spaces) + st.nextToken() + "\n");
		}
		return result;
	}
	
	public static String crop(String string, int length) {
		return string.length() > length ? string.substring(0, length) : string;
	}
	
	public static boolean isFloatNumber(String string) {
		Float v = null;
		try {
			v = Float.parseFloat(string);
		} catch (NumberFormatException e) {
			// ignore
		}		
		return v != null;
	}

	public static boolean isIntegerNumber(String string) {
		Integer v = null;
		try {
			v = Integer.valueOf(string);
		} catch (NumberFormatException e) {
			// ignore
		}		
		return v != null;
	}

	public static boolean containsAnyOf(String toCompare, boolean caseSensitive, String... strings) {
		if (!caseSensitive) {
			toCompare = toCompare.toLowerCase();
		}
		for (String string : strings) {
			if (!caseSensitive) {
				string = string.toLowerCase();
			}
			if (toCompare.contains(string)) {
				return true;
			}
			
		}
		return false;	}
	


	public static boolean isAnyOf(String toCompare, boolean caseSensitive, String... strings) {
		if (!caseSensitive) {
			toCompare = toCompare.toLowerCase();
		}
		for (String string : strings) {
			if (!caseSensitive) {
				string = string.toLowerCase();
			}
			if (toCompare.equals(string)) {
				return true;
			}
			
		}
		return false;
	}

	/**
	 * Creates and returns a String presentation of value with at least digitcount digits (not including negative sign, if any). 
	 */
	public static String toString(int value, int digitCount) {
		String prefix = "";
		String digits;
		
		if (value < 0) {
			// negative
			prefix += "-";
			digits = "" + (-value);
		} else {
			digits = "" + value;
		}
		
		while (digits.length() < digitCount) {
			digits = "0" + digits;
		}
		
		return prefix + digits;
	}
	
	public static String startWithUppercase(String string) {
		if (string.length() > 0) {
			return string.substring(0, 1).toUpperCase() + string.substring(1);
		} else {
			return string;
		}
	}
	
	public static String separateUppercaseChars(String string, String separator) {
		String result = "";
		char prev = 0;
		for (int i = 0; i < string.length(); i++) {
			char current = string.charAt(i);
			if (Character.isUpperCase(current) && Character.isUpperCase(prev)) {
				result += separator;
			} 
			result += current;
			prev = string.charAt(i);
		}
		return result;
	}

}
