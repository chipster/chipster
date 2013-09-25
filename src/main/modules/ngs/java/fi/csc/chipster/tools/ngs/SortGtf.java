package fi.csc.chipster.tools.ngs;

import java.io.File;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.GtfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TsvSorter;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;

public class SortGtf extends JavaAnalysisJobBase {
	
	@Override
	public String getSADL() {
		return 	"TOOL SortGtf.java: \"Sort GTF\" (Sort a GTF file by chromosome and start position.)" + "\n" +
				"INPUT unsorted.gtf: \"GTF file\" TYPE GENERIC" + "\n" +
				"OUTPUT sorted.gtf: \"Sorted GTF file\"" + "\n"; 

	}
	
	
	@Override
	protected void execute() { 
		updateStateToClient(JobState.RUNNING, "sorting");


		try {
			// files
			File inputFile = new File(jobWorkDir, analysis.getInputFiles().get(0).getFileName()); 
			File outputFile = new File(jobWorkDir, analysis.getOutputFiles().get(0).getFileName().getID()); 

			// run sort
			new TsvSorter().sort(
					inputFile, outputFile,
					GtfLineParser.Column.SEQNAME.ordinal(), 
					GtfLineParser.Column.START.ordinal(), new GtfLineParser());

		} catch (Exception e) {
			getResultMessage().setErrorMessage(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED, "");
			return;
		}

		updateStateToClient(JobState.RUNNING, "sorting finished");
	}
}
