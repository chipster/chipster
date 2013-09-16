package fi.csc.chipster.tools.ngs.regions;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.analyser.java.JavaAnalysisHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.RegionOperations;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.util.IOUtils;

public class FindOverlappingDatabaseTool extends RegionTool {

	@Override
	public String getSADL() {
		return 	"TOOL FindOverlappingDatabaseTool.java: \"Match genomic regions against miRBase\" (Returns genomic regions which overlap with miRNA locations reported in miRBase. Adds name of the overlapping miRNA location to the output. Repeating matches are skipped for the TSV output, but included in the BED output.)" + "\n" +
				"INPUT data.bed: \"Set of regions\" TYPE GENERIC" + "\n" +
				"OUTPUT overlapping.bed: \"Overlapping regions\"" + "\n" + 
				"OUTPUT overlapping.tsv: \"Overlapping regions as TSV\"" + "\n" + 
				"PARAMETER database: \"Database\" TYPE [miRBase16: \"miRBase 16\"] DEFAULT miRBase16 (Which miRBase version is used for comparison?)" + 
				"PARAMETER min.overlap.bp: \"Minimum number of overlapping basepairs\" TYPE INTEGER FROM 1 DEFAULT 1 (How many basepairs are required to consider regions overlapping?)";
	}

	@Override
	protected LinkedList<Feature> operate(LinkedList<List<Feature>> inputs, List<String> parameters) throws FileNotFoundException, IOException, URISyntaxException, GBrowserException {

		// Add DB regions to inputs
		File dbDirectory = new File(((JavaAnalysisHandler)this.analysis.getHandler()).getParameters().get("externalToolPath"), "genomic_regions");
		List<Feature> dbRegions = new RegionOperations().loadFile(new File(dbDirectory, "miRBase16.bed"));
		inputs.add(dbRegions);
		
		// Call normal overlapping tool to do actual operation
		LinkedList<String> newParameters = new LinkedList<String>();
		newParameters.add("first_augmented");
		newParameters.add(parameters.get(1)); // min.overlap.bp
		LinkedList<Feature> overlappingRegions = new FindOverlappingTool().operate(inputs, newParameters);
		
		// Do some custom output generation (TSV format)
		FileOutputStream outputStream = null;
		try {
			RegionOperations tool = new RegionOperations();
			outputStream = new FileOutputStream(new File(jobWorkDir, analysis.getOutputFiles().get(1).getFileName().getID())); // overlapping.tsv 
			tool.printTSV(overlappingRegions, outputStream);

		} finally {
			IOUtils.closeIfPossible(outputStream);
		}
		
		// Return overlapping regions for standard output generation 
		return overlappingRegions;
	}
	
	public static void main(String[] args) throws Exception {
		// For testing
		String testPath = "";  
		FindOverlappingDatabaseTool tool = new FindOverlappingDatabaseTool();
		HashMap<String, String> p = new LinkedHashMap<String, String>();
		p.put("externalToolPath", testPath);
		tool.analysis = new ToolDescription(new JavaAnalysisHandler(p));
		tool.analysis.addOutputFile(SADLDescription.Name.createName(testPath + File.separator + "overlapping.bed"), false);
		tool.analysis.addOutputFile(SADLDescription.Name.createName(testPath + File.separator + "overlapping.tsv"), false);
		List<Feature> file1 = new RegionOperations().loadFile(new File(testPath, "miRNAseq.bed"));
		LinkedList<List<Feature>> list = new LinkedList<List<Feature>>();
		list.add(file1);
		LinkedList<String> parameters = new LinkedList<String>();
		parameters.add("");
		parameters.add("1");
		new RegionOperations().print(tool.operate(list, parameters), System.out);
	}
}
