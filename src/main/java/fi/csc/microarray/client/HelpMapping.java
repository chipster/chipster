package fi.csc.microarray.client;

import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.config.DirectoryLayout;

public class HelpMapping {

    public static final String MANUAL_ROOT = DirectoryLayout.getInstance().getConfiguration().getString("client", "manual-root");

	public static String mapToHelppage(OperationDefinition definition) {

		String page = definition.getID();
		if (page.contains(".")) {
			page = page.substring(0, page.lastIndexOf(".")) + ".html"; 
		}
		return MANUAL_ROOT + page;
	}
}
