package fi.csc.chipster.toolbox;

import fi.csc.microarray.description.SADLDescription;

public class ToolboxTool {

	private String id;
	private SADLDescription sadlDescription;
	private String sadlString;
	private String source;
	private String code;

	private String module;
	private String runtime;

	public ToolboxTool() {
	}
	
	public ToolboxTool(String id, SADLDescription sadlDescription, String sadlString, String code, String source, String module, String runtime) {
		this.id = id;
		this.sadlDescription = sadlDescription;
		this.sadlString = sadlString;
		this.source = source;
		this.code = code;
		this.runtime = runtime;
		this.module = module;
	}

	public String getId() {
		return id;
	}
	
	public SADLDescription getSadlDescription() {
		return sadlDescription;
	}
	
	public String getRuntime() {
		return runtime;
	}

	public String getSadlString() {
		return sadlString;
	}

	public String getSource() {
		return source;
	}

	public String getCode() {
		return code;
	}

	public String getModule() {
		return module;
	}	
}
