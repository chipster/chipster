package fi.csc.chipster.toolbox;


public class ToolboxTool {

	private String SADL;
	private String source;
	private String code;

	private String resourceName;
	private String module;
	private String runtime;

	
	public ToolboxTool(String SADL, String code, String source, String resourceName, String module, String runtime) {
		this.SADL = SADL;
		this.source = source;
		this.code = code;
		this.runtime = runtime;
		this.resourceName = resourceName;
		this.module = module;
	}

	public String getRuntime() {
		return runtime;
	}

	public String getResourceName() {
		return resourceName;
	}

	public String getSADL() {
		return SADL;
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
