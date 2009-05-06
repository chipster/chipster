package fi.csc.microarray.analyser;

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
	
	public ParsedRScript parseRScript(InputStream rScriptSource, String commentString) throws MicroarrayException {
		String comment;
		if (commentString == null || commentString.equals("")) {
			comment = COMMENT;
		} else {
			comment = commentString;
		}
		
		ParsedRScript parsedScript = new ParsedRScript();		
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(rScriptSource));
			boolean inHeaderCommentBlock = true;
			for (String line = in.readLine(); line != null; line = in.readLine()) {
				if (!line.startsWith(comment)) {
					inHeaderCommentBlock = false;
				}
				if (inHeaderCommentBlock) {
					parsedScript.VVSADL += line.substring(comment.length());
				}
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

	public ParsedRScript parseRScript(InputStream scriptSource) throws MicroarrayException {
		return parseRScript(scriptSource, COMMENT);
	}

}
