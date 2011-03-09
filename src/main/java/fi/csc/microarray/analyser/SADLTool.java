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

	public static class ParsedScript {
		public String SADL = "";
		public String source = "";
	}
	
	public ParsedScript parseScript(InputStream rScriptSource, String comment) throws IOException {

		ParsedScript parsedScript = new ParsedScript();		
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
				parsedScript.source += line + "\n";
			}
		} finally {
			IOUtils.closeIfPossible(in);
		}
		
		return parsedScript;
	}

}
