package fi.csc.chipster.tools.ngs.regions;

import java.util.LinkedList;
import java.util.List;

import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class FlattenRegionsTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL \"Region operations\" / FlattenRegionsTool.java: \"Flatten overlapping regions\" (Merges overlapping regions of a single file. The returned file does not have any internal overlapping.)" + "\n" +
				"INPUT data1.bed: \"Regions to flatten\" TYPE GENERIC" + "\n" +
				"OUTPUT flattened.bed: \"Merged regions\"" + "\n"; 
	}

	@Override
	protected LinkedList<BpCoordRegion> operate(LinkedList<List<RegionContent>> inputs, List<String> parameters) {
		RegionOperations tool = new RegionOperations();
		return tool.flatten(inputs.get(0));
	}
}
