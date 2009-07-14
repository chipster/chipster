package fi.csc.microarray.util;

public class LangUtil {
	
	public static String simpleClassOf(Object o) {
		try {
			return o.getClass().getSimpleName();
			
		} catch (Exception e) {
			return null;
		}
	}

}
