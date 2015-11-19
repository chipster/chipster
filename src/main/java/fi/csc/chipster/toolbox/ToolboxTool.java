package fi.csc.chipster.toolbox;


public class ToolboxTool {

	private String toolId;
	private String sadl;
	private String source;
	private String code;

	private String resourceName;
	private String module;
	private String runtime;

	public ToolboxTool() {
	}
	
	public ToolboxTool(String toolId, String sadl, String code, String source, String resourceName, String module, String runtime) {
		this.toolId = toolId;
		this.sadl = sadl;
		this.source = source;
		this.code = code;
		this.runtime = runtime;
		this.resourceName = resourceName;
		this.module = module;
	}

	public String getId() {
		return toolId;
	}
	
	public String getRuntime() {
		return runtime;
	}

	public String getResourceName() {
		return resourceName;
	}

	public String getSadl() {
		return sadl;
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
