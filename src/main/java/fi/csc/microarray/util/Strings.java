package fi.csc.microarray.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.apache.commons.lang3.StringUtils;

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

	
	/**
	 * Skip lines from beginning if the line is empty after .trim()
	 * 
	 * Returned string will use "\n" as line break.
	 * @param s
	 * @return
	 * @throws IOException
	 */
	public static String removeEmptyLinesFromBeginning(String s) throws IOException {
		BufferedReader reader = new BufferedReader(new StringReader(s));
		StringBuilder builder = new StringBuilder("");

		try {
			main: for (String line = reader.readLine(); line != null; line = reader.readLine()) {
				
				// if empty line, skip
				if (line.trim().isEmpty()) {
					continue;
				} 
				
				// non-empty line found, get rest of the lines and return
				else {
					builder.append(line + "\n");
					for (String l = reader.readLine(); l != null; l = reader.readLine()) {
						builder.append(l + "\n");
					}
					break main;
				}
			}
		} finally {
			IOUtils.closeIfPossible(reader);
		}

		// remove last line break if original didn't have it
		String result = builder.toString();
		if (result.length() > 0) {
			if (!(s.endsWith("\n") || s.endsWith("\r") || s.endsWith("\r") || s.endsWith("\r\n"))) {
				result = result.substring(0, result.length()-1);
			}
		}
		return result;
	}

	public static String rightPad(String[] array, int width) {
		String line = "";
		for (String title : array) {
			line += StringUtils.rightPad(title, width);
		}
		return line;
	}

	public static List<String> splitConsideringQuotes(String input, char delimiter) {
		// based on http://stackoverflow.com/questions/1757065/splitting-a-comma-separated-string-but-ignoring-commas-in-quotes
		List<String> result = new ArrayList<String>();
		int start = 0;
		boolean inQuotes = false;
		
		for (int current = 0; current < input.length(); current++) {
			
		    if (input.charAt(current) == '\"') {
		    	inQuotes = !inQuotes; // toggle state
		    }
		    
		    boolean atLastChar = (current == input.length() - 1);
		    
		    	
		    if (input.charAt(current) == delimiter && !inQuotes) {
		        result.add(input.substring(start, current).replace("\"", ""));
		        start = current + 1;
		    } else if (atLastChar) {
		    	result.add(input.substring(start).replace("\"", ""));
		    }
		}
		return result;
	}
	
	public static String toHumanReadable(long i) {
		return toHumanReadable(i, true, false);
	}
	
	public static String toRoundedHumanReadable(long i) {
		return toHumanReadable(i, true, true);
	}

	public static String toHumanReadable(long i, boolean returnZero, boolean round) {

		if (i == 0) {
			if (returnZero) {
				return "0";
			} else {
				return "";
			}
		} else if (i < 0) {
			return "" + i;
		}

		int pow = (int) Math.log10(i);

		String sym = "";
		if (pow >= 3) {
			sym = "k";
		}
		if (pow >= 6) {
			sym = "M";
		}
		if (pow >= 9) {
			sym = "G";
		}
		if (pow >= 12) {
			sym = "T";
		}

		long div = (long) Math.pow(10, (pow - pow % 3));		
		
		if (round) {
			// a space between number and symbol to make it less squeezed
			String roundNumber = "" + i / div + " " + sym;
			return roundNumber;
		} else {
					
			String roundNumber = "" + i / div + sym;
			String remainder = toHumanReadable(i % div, false, false);
			return  roundNumber + " " + remainder;
		}
	}

	/**
	 * Split a string using a method String.split() if the string is not empty.
	 * In case of empty string, return a zero length array as opposed to
	 * String.split(), which returns an array of length one containing an empty
	 * string.
	 * 
	 * @param str
	 * @param regex
	 * @return
	 */
	public static String[] splitUnlessEmpty(String str,
			String regex) {
		if (str == null || str.isEmpty()) {
			return new String[0];
		} else {
			return str.split(regex);
		}
	}
}
