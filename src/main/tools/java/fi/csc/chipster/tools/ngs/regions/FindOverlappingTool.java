package fi.csc.chipster.tools.ngs.regions;

import java.util.LinkedList;
import java.util.List;

import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class FindOverlappingTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL \"Region operations\" / FindOverlappingTool.java: \"Find overlapping regions\" (Returns regions that have overlap with some region in the other input file. Also known as intersection.)" + "\n" +
				"INPUT data1.bed: \"First set of regions\" TYPE GENERIC" + "\n" +
				"INPUT data2.bed: \"Second set of regions\" TYPE GENERIC" + "\n" +
				"OUTPUT overlapping.bed: \"Overlapping regions\"" + "\n" + 
				"PARAMETER return.type: \"Type of returned regions\" TYPE [merged: \"Merged overlapping regions\", intersection: \"Overlapping piece of regions\", original: \"Both of the original regions\"] DEFAULT merged (How overlapping regions are returned?)" + 
				"PARAMETER min.overlap.bp: \"Minimum number of overlapping basepairs\" TYPE INTEGER FROM 1 DEFAULT 1 (How many basepairs are required to consider regions overlapping?)";
	}

	@Override
	protected LinkedList<BpCoordRegion> operate(LinkedList<List<RegionContent>> inputs, List<String> parameters) {
		RegionOperations tool = new RegionOperations();
		RegionOperations.PairPolicy pairPolicy;
		if ("intersection".equals(inputMessage.getParameters().get(0))) {
			pairPolicy = RegionOperations.INTERSECT_PAIR_POLICY; 
		} else if ("original".equals(inputMessage.getParameters().get(0))) {
			pairPolicy = RegionOperations.ORIGINALS_PAIR_POLICY;
		} else {
			pairPolicy = RegionOperations.MERGE_PAIR_POLICY;
		}
		Long minOverlap = Long.parseLong(inputMessage.getParameters().get(1));
		return tool.intersect(inputs.get(0), inputs.get(1), minOverlap, pairPolicy);
		
	}
}
