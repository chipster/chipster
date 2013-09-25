package fi.csc.chipster.tools.ngs.regions;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.RegionOperations;

public class CombineRegionsTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL CombineRegionsTool.java: \"Combine region files\" (Returns combined regions from both input files. Also known as union.)" + "\n" +
				"INPUT data1.bed: \"Region file A\" TYPE GENERIC" + "\n" +
				"INPUT data2.bed: \"Region file B\" TYPE GENERIC" + "\n" +
				"OUTPUT combined.bed: \"Combined regions\"" + "\n" + 
				"PARAMETER merge.overlapping: \"Merge overlapping regions before returning them\" TYPE [yes: \"Yes\", no: \"No\"] DEFAULT yes (Should result be flattened?)" +
				"PARAMETER min.overlap.bp: \"Minimum number of overlapping bases, if merging\" TYPE INTEGER FROM 1 DEFAULT 1 (If result is flattened, how many bases are required to consider regions overlapping?)";
	}

	@Override
	protected LinkedList<Feature> operate(LinkedList<List<Feature>> inputs, List<String> parameters) {
		RegionOperations tool = new RegionOperations();
		boolean flatten = "yes".equals(parameters.get(0));
		Long minOverlap = Long.parseLong(parameters.get(1));
		return tool.merge(inputs.get(0), inputs.get(1), minOverlap, flatten);
	}
}
