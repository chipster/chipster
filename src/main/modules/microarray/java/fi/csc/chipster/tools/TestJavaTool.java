package fi.csc.chipster.tools;

import java.io.File;
import java.io.FileWriter;

import org.apache.commons.io.FileUtils;

import fi.csc.microarray.comp.JobCancelledException;
import fi.csc.microarray.comp.java.JavaCompJobBase;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;

public class TestJavaTool extends JavaCompJobBase {

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
			commentWriter.write(inputMessage.getParameters(JAVA_PARAMETER_SECURITY_POLICY, this.toolDescription).get(0));
			commentWriter.flush();
			commentWriter.close();

		} catch (Exception ioe) {
			outputMessage.setErrorMessage("Running Java job failed.");
			outputMessage.setOutputText(Exceptions.getStackTrace(ioe));
			updateState(JobState.FAILED, "");
		}
		updateState(JobState.RUNNING, "Java tool finished");
	}

}
