package fi.csc.chipster.tools.ngs.regions;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.RegionOperations;

public class FuseRegionsTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL FuseRegionsTool.java: \"Fuse overlapping regions\" (Merges overlapping regions of a single file. The returned file does not have any internal overlapping.)" + "\n" +
				"INPUT data.bed: \"Regions to fuse\" TYPE GENERIC" + "\n" +
				"OUTPUT fused.bed: \"Merged regions\"" + "\n"; 
	}

	@Override
	protected LinkedList<Feature> operate(LinkedList<List<Feature>> inputs, List<String> parameters) {
		RegionOperations tool = new RegionOperations();
		return tool.flatten(inputs.get(0));
	}
}
