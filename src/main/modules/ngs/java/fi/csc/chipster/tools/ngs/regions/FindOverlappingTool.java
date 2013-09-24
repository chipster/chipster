package fi.csc.chipster.tools.ngs.regions;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.RegionOperations;

public class FindOverlappingTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL FindOverlappingTool.java: \"Find overlapping regions\" (Returns regions that have overlap with some region in the other input file. Also known as intersection.)" + "\n" +
				"INPUT data1.bed: \"Region file A\" TYPE GENERIC" + "\n" +
				"INPUT data2.bed: \"Region file B\" TYPE GENERIC" + "\n" +
				"OUTPUT overlapping.bed: \"Overlapping regions\"" + "\n" + 
				"PARAMETER return.type: \"Type of returned regions\" TYPE [first: \"Original regions from file A\", first_augmented: \"Original regions from file A augmented with info from file B\", both: \"Original regions from both files\", merged: \"Regions from both files merged\", intersection: \"Overlapping parts of the regions\"] DEFAULT first (How overlapping regions are returned?)" + 
				"PARAMETER min.overlap.bp: \"Minimum number of overlapping bases\" TYPE INTEGER FROM 1 DEFAULT 1 (How many bases are required to consider regions overlapping?)";
	}
	

	@Override
	protected LinkedList<Feature> operate(LinkedList<List<Feature>> inputs, List<String> parameters) {
		RegionOperations tool = new RegionOperations();
		RegionOperations.PairPolicy pairPolicy;
		if ("intersection".equals(parameters.get(0))) {
			pairPolicy = RegionOperations.INTERSECT_PAIR_POLICY; 
		} else if ("both".equals(parameters.get(0))) {
			pairPolicy = RegionOperations.ORIGINALS_PAIR_POLICY;
		} else if ("merged".equals(parameters.get(0))) {
			pairPolicy = RegionOperations.MERGE_PAIR_POLICY;
		} else if ("first_augmented".equals(parameters.get(0))) {
			pairPolicy = RegionOperations.LEFT_PAIR_POLICY_WITH_AUGMENTATION;
		} else {
			pairPolicy = RegionOperations.LEFT_PAIR_POLICY;
		}
		Long minOverlap = Long.parseLong(parameters.get(1));
		return tool.intersect(inputs.get(0), inputs.get(1), minOverlap, pairPolicy, false);
		
	}
}
