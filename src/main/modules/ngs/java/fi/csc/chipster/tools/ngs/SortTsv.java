package fi.csc.chipster.tools.ngs;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.TsvLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TsvSorter;
import fi.csc.microarray.comp.java.JavaCompJobBase;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;

public class SortTsv extends JavaCompJobBase {

	public static final String COLUMN_ID = "column";
	public static final String FIRST_ID = "first";
	public static final String SECOND_ID = "second";
	
	@Override
	public String getSADL() {
		return 	"TOOL SortTsv.java: \"Sort TSV\" (Sort a TSV file by chromosome and start position.)" + "\n" +
				"INPUT input.tsv: \"TSV file\" TYPE GENERIC" + "\n" +
				"OUTPUT sorted.tsv: \"Sorted TSV file\"" + "\n" + 
				"PARAMETER " + COLUMN_ID + ": \"Chromosome column\" TYPE [" + FIRST_ID + ": First, " + SECOND_ID + ": Second] DEFAULT " + FIRST_ID + " (Select the column that contains chromosome information.)" + "\n"; 
	}
	
	
	@Override
	protected void execute() { 
		updateStateToClient(JobState.RUNNING, "sorting");

		try {
			// files
			File inputFile = new File(jobDataDir, toolDescription.getInputFiles().get(0).getFileName()); 
			File outputFile = new File(jobDataDir, toolDescription.getOutputFiles().get(0).getFileName().getID());
			
			List<String> parameters = inputMessage.getParameters(JAVA_PARAMETER_SECURITY_POLICY, toolDescription);			
			String columnString = parameters.get(0);

			int chrColumn = 0; //if not second, this will apply
			if (SECOND_ID.equals(columnString)) {
				chrColumn = 1;
			}
			
			sort(inputFile, outputFile, chrColumn);	

		} catch (Exception e) {			
			
			getResultMessage().setErrorMessage(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED, "");
			return;
		}

		updateStateToClient(JobState.RUNNING, "sorting finished");
	}
	
	private static void sort(File inputFile, File outputFile, int chrColumn) throws MalformedURLException, IOException, URISyntaxException, Exception {
		// run sort
		new TsvSorter().sort(
				inputFile, outputFile, 
				chrColumn, chrColumn + 1, new TsvLineParser(new DataUrl(inputFile), chrColumn));
	}


	public static void main(String[] args) throws Exception {

		try {

			File in = new File(args[0]);
			File out = new File(args[1]);
			int chrColumn = Integer.parseInt(args[2]) - 1;

			sort(in, out, chrColumn);

		} catch (Exception e) {
			e.printStackTrace();
						
			System.out.println(
					"usage: SortGtf input-file output-file chr-column\n" +
					"  chr-column:	number of column containin chromosomes, starting from 1\n" +
					"  \n" + 
					"example:\n " +
					"  java -cp chipster-3.0.0.jar fi.csc.chipster.tools.ngs.SortTsv de-genes-cufflinks.tsv de-genes-cufflinks-sort.tsv 2");
		}				
	}
}
