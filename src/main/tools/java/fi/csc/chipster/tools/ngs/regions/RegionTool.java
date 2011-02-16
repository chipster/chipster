package fi.csc.chipster.tools.ngs.regions;

import java.io.File;
import java.io.FileOutputStream;
import java.util.LinkedList;
import java.util.List;

import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.IOUtils;

public abstract class RegionTool extends JavaAnalysisJobBase {

	protected abstract LinkedList<RegionContent> operate(LinkedList<List<RegionContent>> inputs, List<String> parameters);
	
	@Override
	protected void execute() { 
		try {
			updateStateToClient(JobState.RUNNING, "preprocessing");

			// Parse inputs
			RegionOperations tool = new RegionOperations();
			LinkedList<List<RegionContent>> inputs = new LinkedList<List<RegionContent>>();
			for (int i = 0; i < analysis.getInputFiles().size(); i++) {
				File inputFile = new File(jobWorkDir, analysis.getInputFiles().get(i).getFileName());
				inputs.add(tool.loadFile(inputFile));
			}

			// Delegate actual processing to subclasses
			LinkedList<RegionContent> output = operate(inputs, inputMessage.getParameters());
			
			// Sort result
			new RegionOperations().sort(output);
			
			// Write output
			FileOutputStream outputStream = null;
			try {
				outputStream = new FileOutputStream(new File(jobWorkDir, analysis.getOutputFiles().get(0).getFileName().getID())); 
				tool.printRegions(output, outputStream);

			} finally {
				IOUtils.closeIfPossible(outputStream);
			}
			
		} catch (Exception e) {
			updateState(JobState.FAILED, e.getMessage());
			return;
		}
		updateStateToClient(JobState.RUNNING, "preprocessing finished");
	}
	


}
