package fi.csc.chipster.toolbox.toolpartsparser;

import java.io.File;

import fi.csc.chipster.toolbox.SADLTool.ParsedScript;

public interface ToolPartsParser {

	public ParsedScript parse(File moduleDir, String toolFilename) throws Exception;

}
