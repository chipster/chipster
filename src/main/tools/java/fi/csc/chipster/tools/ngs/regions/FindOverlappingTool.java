package fi.csc.chipster.tools.ngs.regions;

import java.util.LinkedList;
import java.util.List;

import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class FindOverlappingTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL \"Region operations\" / FindOverlappingTool.java: \"Find overlapping regions\" (Returns regions that have overlap with some region in the other input file. Also known as intersection.)" + "\n" +
				"INPUT data1.bed: \"First set of regions\" TYPE GENERIC" + "\n" +
				"INPUT data2.bed: \"Second set of regions\" TYPE GENERIC" + "\n" +
				"OUTPUT overlapping.bed: \"Overlapping regions\"" + "\n" + 
				"PARAMETER return.type: \"Type of returned regions\" TYPE [first: \"Regions from the first set only\", both: \"Regions from both sets\", merged: \"Merged regions\", intersection: \"Overlapping pieces of regions\"] DEFAULT first (How overlapping regions are returned?)" + 
				"PARAMETER min.overlap.bp: \"Minimum number of overlapping basepairs\" TYPE INTEGER FROM 1 DEFAULT 1 (How many basepairs are required to consider regions overlapping?)";
	}

	@Override
	protected LinkedList<RegionContent> operate(LinkedList<List<RegionContent>> inputs, List<String> parameters) {
		RegionOperations tool = new RegionOperations();
		RegionOperations.PairPolicy pairPolicy;
		if ("intersection".equals(parameters.get(0))) {
			pairPolicy = RegionOperations.INTERSECT_PAIR_POLICY; 
		} else if ("both".equals(parameters.get(0))) {
			pairPolicy = RegionOperations.ORIGINALS_PAIR_POLICY;
		} else if ("merged".equals(parameters.get(0))) {
			pairPolicy = RegionOperations.MERGE_PAIR_POLICY;
		} else {
			pairPolicy = RegionOperations.LEFT_PAIR_POLICY;
		}
		Long minOverlap = Long.parseLong(parameters.get(1));
		return tool.intersect(inputs.get(0), inputs.get(1), minOverlap, pairPolicy, false);
		
	}
}
