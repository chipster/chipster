package fi.csc.chipster.tools.ngs.regions;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class CombineRegionsTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL \"Region operations\" / CombineRegionsTool.java: \"Combine region files\" (Returns combined regions from both of the inputs. Also known as union.)" + "\n" +
				"INPUT data1.bed: \"First input file\" TYPE GENERIC" + "\n" +
				"INPUT data2.bed: \"Second input file\" TYPE GENERIC" + "\n" +
				"OUTPUT result.bed: \"Combined regions\"" + "\n" + 
				"PARAMETER merge.overlapping: \"Merge overlapping regions before returning them\" TYPE [yes: \"Yes\", no: \"No\"] DEFAULT yes (Should result be flattened?)" +
				"PARAMETER min.overlap.bp: \"Minimum number of overlapping basepairs, if merging\" TYPE INTEGER FROM 1 DEFAULT 1 (If result is flattened, how many basepairs are required to consider regions overlapping?)";
	}

	@Override
	protected LinkedList<BpCoordRegion> operate(LinkedList<List<RegionContent>> inputs, List<String> parameters) {
		RegionOperations tool = new RegionOperations();
		boolean flatten = "yes".equals(inputMessage.getParameters().get(0));
		Long minOverlap = Long.parseLong(inputMessage.getParameters().get(1));
		return tool.merge(inputs.get(0), inputs.get(1), minOverlap, flatten);
	}
}
