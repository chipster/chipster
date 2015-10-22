package fi.csc.microarray.comp;

import fi.csc.chipster.toolbox.Toolbox;
import fi.csc.chipster.toolbox.ToolboxClient;
import fi.csc.chipster.toolbox.ToolboxTool;

public class OldToolboxClient implements ToolboxClient {

	private Toolbox toolbox;
	
	public OldToolboxClient(Toolbox toolbox) {
		this.toolbox = toolbox;
	}
	
	@Override
	public ToolboxTool getTool(String id) {
		return toolbox.getTool(id);
	}

}
