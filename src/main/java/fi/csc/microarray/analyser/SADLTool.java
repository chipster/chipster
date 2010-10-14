package fi.csc.microarray.analyser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import fi.csc.microarray.util.IOUtils;

/**
 * Tool for SADL related tasks.
 * 
 * @author Aleksi Kallio
 */
public class SADLTool {

	/**
	 * The String that is used to begin comment lines in scripts. For parsing
	 * SADL from comments.
	 */
	private static String COMMENT = "#";

	public static class ParsedRScript {
		public String SADL = "";
		public String rSource = "";
	}
	
	public ParsedRScript parseRScript(InputStream rScriptSource, String commentString) throws IOException {
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
					String strippedLine = line.substring(comment.length());
					if (strippedLine.startsWith(" ")) {
						strippedLine = strippedLine.substring(1); // remove space after comment symbol
					}
					parsedScript.SADL += strippedLine + "\n";
				}
				parsedScript.rSource += line + "\n";
			}
		} finally {
			IOUtils.closeIfPossible(in);
		}
		return parsedScript;
	}

	public ParsedRScript parseRScript(InputStream scriptSource) throws IOException {
		return parseRScript(scriptSource, COMMENT);
	}

}
