package fi.csc.chipster.toolbox;

import fi.csc.chipster.toolbox.SADLTool.ParsedScript;

public class ToolboxTool {

	private ParsedScript parsedScript;
	private String resourceName;
	private String module;
	private String runtime;

	
	public ToolboxTool(ParsedScript parsedScript, String resourceName, String module, String runtime) {
		this.parsedScript = parsedScript;
		this.runtime = runtime;
		this.resourceName = resourceName;
		this.module = module;
	}

	public ParsedScript getParsedScript() {
		return parsedScript;
	}

	public String getRuntime() {
		return runtime;
	}

	public String getResourceName() {
		return resourceName;
	}

	public String getModule() {
		return module;
	}

	
	
}
