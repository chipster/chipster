package fi.csc.chipster.tools.ngs;

import java.io.File;
import java.util.List;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.TsvLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TsvSorter;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;

public class SortTsv extends JavaAnalysisJobBase {

	public static final String COLUMN_ID = "column";
	public static final String FIRST_ID = "first";
	public static final String SECOND_ID = "second";
	
	@Override
	public String getSADL() {
		return 	"TOOL SortTsv.java: \"Sort TSV\" (Sort a TSV file by chromosome and start position.)" + "\n" +
				"INPUT input.tsv: \"TSV file\" TYPE GENERIC" + "\n" +
				"OUTPUT sorted.tsv: \"Sorted TSV file\"" + "\n" + 
				"PARAMETER " + COLUMN_ID + ": \"Chromosome column\" TYPE [" + FIRST_ID + ": First, " + SECOND_ID + ": Second] DEFAULT First (Select the column that contains chromosome information.)" + "\n"; 
	}
	
	
	@Override
	protected void execute() { 
		updateStateToClient(JobState.RUNNING, "sorting");

		try {
			// files
			File inputFile = new File(jobWorkDir, analysis.getInputFiles().get(0).getFileName()); 
			File outputFile = new File(jobWorkDir, analysis.getOutputFiles().get(0).getFileName().getID());
			
			List<String> parameters = inputMessage.getParameters(JAVA_PARAMETER_SECURITY_POLICY, analysis);			
			String columnString = parameters.get(0);

			int chrColumn = 0; //if not second, this will apply
			if (SECOND_ID.equals(columnString)) {
				chrColumn = 1;
			}
			
			// run sort
			new TsvSorter().sort(
					inputFile, outputFile, 
					chrColumn, chrColumn + 1, new TsvLineParser(new DataUrl(inputFile), chrColumn));

		} catch (Exception e) {			
			
			getResultMessage().setErrorMessage(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED, "");
			return;
		}

		updateStateToClient(JobState.RUNNING, "sorting finished");
	}
}
