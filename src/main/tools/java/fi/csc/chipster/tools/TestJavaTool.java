package fi.csc.chipster.tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.io.FileUtils;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.messaging.JobState;

public class TestJavaTool extends JavaAnalysisJobBase {

	@Override
	public String getSADL() {
		return 	" ANALYSIS Test/JavaTool (Simple JavaTool test.)" + "\n" + 
				" INPUT GENERIC input.tsv OUTPUT output.tsv, comment.txt" + "\n" +
				" PARAMETER hyvinkö.menee [yes, no] DEFAULT yes (Hyvä parametri)";
 
	}

	@Override
	protected void execute() throws JobCancelledException {
		updateStateDetailToClient("Java tool running");
		
		File inputFile = new File(jobWorkDir, "input.tsv");
		File outputFile = new File(jobWorkDir, "output.tsv");
		try {
			FileUtils.copyFile(inputFile, outputFile);

			File commentFile = new File(jobWorkDir, "comment.txt");
			FileWriter commentWriter = new FileWriter(commentFile);
			commentWriter.write(inputMessage.getParameters().get(0));
			commentWriter.flush();
			commentWriter.close();

		} catch (IOException ioe) {
			outputMessage.setErrorMessage("Running Java job failed.");
			outputMessage.setOutputText(ioe.toString());
			updateState(JobState.FAILED, "");
		}
		updateState(JobState.RUNNING, "Java tool finished");
	}

}
