package fi.csc.microarray.util;

import java.util.StringTokenizer;

public class Strings {

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

}
