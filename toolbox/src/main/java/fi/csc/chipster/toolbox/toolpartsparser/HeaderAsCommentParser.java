package fi.csc.chipster.toolbox.toolpartsparser;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;

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
	public ParsedScript parse(Path moduleDir, String toolFilename) throws IOException {

		Path toolFile = moduleDir.resolve(toolPath).resolve(toolFilename);

		SADLTool.ParsedScript parsedScript;
		try (InputStream scriptSource = Files.newInputStream(toolFile)) {

			// read the SADL from the comment block in the beginning of file
			// and the actual source code
			parsedScript = new SADLTool(commentString).parseScript(scriptSource);
		}
		return parsedScript;
	}
}
