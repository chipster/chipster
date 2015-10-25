package fi.csc.microarray.comp;

import java.util.LinkedList;

import org.junit.Before;

import fi.csc.microarray.client.RemoteServiceAccessor;
import fi.csc.microarray.client.ServiceAccessor;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class ToolTestBase extends AnalysisTestBase {
	
	protected ServiceAccessor serviceAccessor;
	protected LinkedList<ToolModule> toolModules = new LinkedList<ToolModule>();
	
	public ToolTestBase(String username, String password, String configURL) {
		super(username, password, configURL);
	}

	@Before
	public void setUp() throws Exception {
		super.setUp();
		serviceAccessor = new RemoteServiceAccessor();
		serviceAccessor.initialise(this.manager, this.authenticationListener);
		serviceAccessor.fetchDescriptions(new MicroarrayModule());
		toolModules.addAll(serviceAccessor.getModules());
	}

	/**
	 * 
	 * @param toolId
	 * @return null if operation definition is not found
	 */
	protected OperationDefinition getOperationDefinition(String toolId) {
		for (ToolModule module : toolModules) {
			OperationDefinition tool = module.getOperationDefinition(toolId);
			if (tool != null) {
				return tool;
			}
		}
		return null;
	}
}



