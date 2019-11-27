package fi.csc.chipster.toolbox.toolpartsparser;

import java.nio.file.Path;

import fi.csc.chipster.toolbox.SADLTool.ParsedScript;

public interface ToolPartsParser {

	public ParsedScript parse(Path moduleDir, String toolFilename) throws Exception;

}
