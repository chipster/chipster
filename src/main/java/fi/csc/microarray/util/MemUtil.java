package fi.csc.microarray.util;

import java.io.PrintWriter;
import java.io.StringWriter;

public class MemUtil {

	public static String bytesToMegas(long bytes) {
		float megas = ((float)bytes)/1024f/1024f;
		StringWriter s = new StringWriter(); 
		new PrintWriter(s).printf("%.0f", megas);
		return s.toString();
	}
	
	public static long getUsed() {
		return Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
	}
	
	public static String getMemInfo() {
		long used = getUsed();
		//String totalS = bytesToMegas(Runtime.getRuntime().totalMemory());
		String usedS = bytesToMegas(used);
		String maxS = bytesToMegas(Runtime.getRuntime().maxMemory());
		
		return usedS + "M / " + maxS + "M";
	}

}
