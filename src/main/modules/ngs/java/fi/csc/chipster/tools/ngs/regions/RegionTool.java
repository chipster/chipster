package fi.csc.chipster.tools.ngs.regions;

import java.io.File;
import java.io.FileOutputStream;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.RegionOperations;
import fi.csc.microarray.comp.java.JavaCompJobBase;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.IOUtils;

public abstract class RegionTool extends JavaCompJobBase {

	protected abstract LinkedList<Feature> operate(LinkedList<List<Feature>> inputs, List<String> parameters) throws Exception;
	
	@Override
	protected void execute() { 
		try {
			updateStateToClient(JobState.RUNNING, "preprocessing");

			// Parse inputs
			RegionOperations tool = new RegionOperations();
			LinkedList<List<Feature>> inputs = new LinkedList<List<Feature>>();
			for (int i = 0; i < toolDescription.getInputFiles().size(); i++) {
				File inputFile = new File(jobWorkDir, toolDescription.getInputFiles().get(i).getFileName());
				inputs.add(tool.loadFile(inputFile));
			}

			// Delegate actual processing to subclasses
			List<String> parameters = inputMessage.getParameters(JAVA_PARAMETER_SECURITY_POLICY, toolDescription);
			LinkedList<Feature> output = operate(inputs, parameters);
			
			// Sort result
			new RegionOperations().sort(output);
			
			// Write output
			FileOutputStream outputStream = null;
			try {
				outputStream = new FileOutputStream(new File(jobWorkDir, toolDescription.getOutputFiles().get(0).getFileName().getID())); 
				tool.print(output, outputStream);

			} finally {
				IOUtils.closeIfPossible(outputStream);
			}
			
		} catch (Exception e) {
			updateState(JobState.FAILED, e.getMessage());
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
			return;
		}
		updateStateToClient(JobState.RUNNING, "preprocessing finished");
	}
	


}
