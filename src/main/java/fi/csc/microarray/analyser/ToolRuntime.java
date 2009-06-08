package fi.csc.microarray.analyser;

public class ToolRuntime {

	private String name;
	private AnalysisHandler handler;
	private boolean disabled;

	public ToolRuntime(String name, AnalysisHandler handler, boolean disabled) {
		this.name = name;
		this.handler = handler;
		this.disabled = disabled;
	}

	public String getName() {
		return name;
	}

	public AnalysisHandler getHandler() {
		return handler;
	}

	public boolean isDisabled() {
		return disabled || handler.isDisabled();
	}

}
