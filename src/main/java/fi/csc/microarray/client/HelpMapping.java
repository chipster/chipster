package fi.csc.microarray.client;

import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.config.DirectoryLayout;

public class HelpMapping {

    public static final String MANUAL_ROOT = DirectoryLayout.getInstance().getConfiguration().getString("client", "manual-root");

	public static String mapToHelppage(OperationDefinition definition) {
		
		// blast
		if ("BLAST".equals(definition.getCategory().getName())) {
			return "http://www.csc.fi/english/research/sciences/bioscience/programs/blast/index_html"; 
		}
		
		// mafft
		else if ("Alignment:Multiple".equals(definition.getCategory().getName()) && definition.getID().startsWith("mafft")) {
			return "http://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html";

		}
		
		// others
		else {
			String page = definition.getID();
			if (page.contains(".")) {
				page = page.substring(0, page.lastIndexOf(".")) + ".html"; 
			}
			return MANUAL_ROOT + page;
		}
	}
}
