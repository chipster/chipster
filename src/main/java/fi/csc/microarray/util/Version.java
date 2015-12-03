package fi.csc.microarray.util;

import java.util.LinkedList;
import java.util.List;

/**
 * Class for storing and comparing version number.
 * 
 * Format: xx.yyy.z.hhhh
 * 
 * Parts separated by '.', number of parts and their length not limited.
 * 
 */
public class Version {

	private Integer[] parts;

	private Version(Integer[] parts) {
		this.parts = parts;
	}

    public static Version parse(String version) {
    	String[] strings = version.split("\\.", -1);

    	List<Integer> list = new LinkedList<Integer>();
    	for (String s: strings) {
    		list.add(Integer.parseInt(s));
    	}
    	
        return new Version(list.toArray(new Integer[] {}));
    }

	
	public int compareTo(Version v) {
		int minLength = parts.length < v.getParts().length ? parts.length : v.getParts().length;
		
		for (int i = 0; i < minLength; i++) {
			if (parts[i] != v.getParts()[i]) {
				return parts[i].compareTo(v.getParts()[i]);
			}
		}
		
		// gone through all the parts of the shorter version, all being same
		// return longer, or 0 if the same length
		return Integer.compare(parts.length, v.getParts().length);
	}


	Integer[] getParts() {
		return this.parts;
	}


	@Override
	public String toString() {
		String s = "";
		for (Integer i: parts) {
			s += i + ".";
		}
		if (s.endsWith(".")) {
			s = s.substring(0, s.length()-1);
		}

		return s;
	}
}
