package fi.csc.chipster.toolbox.toolpartsparser;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;

import de.schlichtherle.truezip.file.TFile;
import de.schlichtherle.truezip.file.TFileInputStream;
import fi.csc.chipster.toolbox.SADLTool;
import fi.csc.chipster.toolbox.SADLTool.ParsedScript;

public class HeaderAsCommentParser implements ToolPartsParser {

	private String commentString;
	private String toolPath;
	
	public HeaderAsCommentParser(String commentString, String toolPath) {
		this.commentString = commentString;
		this.toolPath = toolPath;
	}
	
	@Override
	public ParsedScript parse(File moduleDir, String toolFilename) throws IOException {

		TFile toolFile = new TFile(moduleDir, toolPath + File.separator + toolFilename);
		InputStream scriptSource;
		scriptSource = new TFileInputStream(toolFile);

		// read the SADL from the comment block in the beginning of file
		// and the actual source code
		SADLTool.ParsedScript parsedScript;
		parsedScript = new SADLTool(commentString).parseScript(scriptSource);

		return parsedScript;
	}

}
