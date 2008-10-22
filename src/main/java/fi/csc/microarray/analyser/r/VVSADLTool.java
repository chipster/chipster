package fi.csc.microarray.analyser.r;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import fi.csc.microarray.MicroarrayException;

/**
 * Tool for VVSADL related tasks.
 * 
 * @author Aleksi Kallio
 */
public class VVSADLTool {

	/**
	 * The String that is used to begin comment lines in scripts. For parsing
	 * VVSADL from comments.
	 */
	private static String COMMENT = "#";

	public static class ParsedRScript {
		public String VVSADL = "";
		public String rSource = "";
	}
	
	public ParsedRScript parseRScript(InputStream rScriptSource) throws MicroarrayException {
		ParsedRScript parsedScript = new ParsedRScript();		
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(rScriptSource));
			for (String line = in.readLine(); line != null && line.startsWith(COMMENT); line = in.readLine()) {
				parsedScript.VVSADL += line.substring(COMMENT.length());
			}
			for (String line = in.readLine(); line != null; line = in.readLine()) {
				parsedScript.rSource += line + "\n";
			}
		} catch (IOException e) {
			throw new MicroarrayException(e);
		} finally {
			try {
				in.close();
			} catch (Exception e) {
				// ignore
			}
		}
		return parsedScript;
	}
}
