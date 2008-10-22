package fi.csc.microarray.util;

public class Exceptions {

	public static String getStackTrace(Throwable throwable) {
	
		String trace = throwable.getClass().getSimpleName();
		if (throwable.getMessage() != null) {
			trace += ": " + throwable.getMessage();
		}
		trace += "\n";

		for (StackTraceElement element : throwable.getStackTrace()) {
			trace += (element.toString() + "\n"); 
		}
		return trace;
	}
}
