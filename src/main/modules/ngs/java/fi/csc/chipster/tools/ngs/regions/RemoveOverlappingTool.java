package fi.csc.chipster.tools.ngs.regions;

import java.util.LinkedList;
import java.util.List;

import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class RemoveOverlappingTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL RemoveOverlappingTool.java: \"Remove overlapping regions\" (Returns regions in the first input file that do not have overlap with any of the regions in the second input file. Also known as subtraction.)" + "\n" +
				"INPUT data1.bed: \"Regions to remove from\" TYPE GENERIC" + "\n" +
				"INPUT data2.bed: \"Regions that possibly overlap with the first input file\" TYPE GENERIC" + "\n" +
				"OUTPUT nonoverlapping.bed: \"Regions of first input that do not overlap\"" + "\n" + 
				"PARAMETER min.overlap.bp: \"Minimum number of overlapping basepairs\" TYPE INTEGER FROM 1 DEFAULT 1 (How many basepairs are required to consider regions overlapping?)";
	}

	@Override
	protected LinkedList<RegionContent> operate(LinkedList<List<RegionContent>> inputs, List<String> parameters) {
		RegionOperations tool = new RegionOperations();
		Long minOverlap = Long.parseLong(parameters.get(0));
		return tool.subtract(inputs.get(0), inputs.get(1), minOverlap);
	}
}
