package fi.csc.microarray.util;

import java.io.IOException;

import javax.net.ssl.SSLHandshakeException;

import fi.csc.microarray.exception.MicroarrayException;


public class Exceptions {

	public static String getStackTrace(Throwable throwable) {
	
		String trace = throwable.getClass().getSimpleName();
		if (throwable.getMessage() != null) {
			trace += ": " + throwable.getMessage();
		}
		trace += "\n";

		trace += getStackTrace(throwable.getStackTrace());
		
		if (throwable.getCause() != null) {
			trace += "Caused by: ";
			trace += getStackTrace(throwable.getCause());
		}
		return trace;
	}
	
	public static String getStackTrace(StackTraceElement[] stackTrace) {
		String trace = "";
		for (StackTraceElement element : stackTrace) {
			trace += (element.toString() + "\n"); 
		}
		return trace;
	}

	/**
	 * Check if this exception or any of its causes is a instance of the given
	 * type.
	 * 
	 * Type must be given as an object (i.e. new IOException("")), because
	 * that's what the class.isInstance() expects.
	 * 
	 * @param exception
	 * @param type
	 * @return
	 */
	public static boolean isCausedBy(Throwable exception, Object type) {
		
		if (type.getClass().isInstance(exception)) {
			return true;
		}
		
		Throwable e = exception;
		
		while (e.getCause() != null) {
			e = e.getCause();
			
			if (type.getClass().isInstance(e)) {
				return true;
			}
		}
		return false;
	}
	
	public static void main(String args[]) throws InstantiationException, IllegalAccessException {
		
		System.out.println(isCausedBy(new MicroarrayException(new Exception()), new SSLHandshakeException("")));
		System.out.println(isCausedBy(new MicroarrayException(new SSLHandshakeException("")), new SSLHandshakeException("")));
		
		System.out.println(isCausedBy(new MicroarrayException(""), new SSLHandshakeException("")));
		System.out.println(isCausedBy(new SSLHandshakeException(""), new IOException("")));			
	}
}
