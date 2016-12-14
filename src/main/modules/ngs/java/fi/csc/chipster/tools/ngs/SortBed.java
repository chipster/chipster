package fi.csc.chipster.tools.ngs;

import java.io.File;

import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.BedLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TsvSorter;
import fi.csc.microarray.comp.java.JavaCompJobBase;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;

public class SortBed extends JavaCompJobBase {
	
	@Override
	public String getSADL() {
		return 	"TOOL fi.csc.chipster.tools.ngs.SortBed.java: \"Sort BED\" (Sort a BED file by chromosome and start position.)" + "\n" +
				"INPUT regions.bed: \"BED file\" TYPE GENERIC" + "\n" +
				"OUTPUT sorted.bed: \"Sorted BED file\"" + "\n"; 

	}
	
	
	@Override
	protected void execute() { 
		updateStateToClient(JobState.RUNNING, "sorting");


		try {
			// files
			File inputFile = new File(jobDataDir, toolDescription.getInputFiles().get(0).getFileName()); 
			File outputFile = new File(jobDataDir, toolDescription.getOutputFiles().get(0).getFileName().getID()); 

			// run sort
			//BEDParser increments coordinates by one, but it's not a problem because only its column order is used
			new TsvSorter().sort(
					inputFile, outputFile,
					BedLineParser.Column.CHROM.ordinal(), BedLineParser.Column.CHROM_START.ordinal(), new BedLineParser(false));

		} catch (Exception e) {
			getResultMessage().setErrorMessage(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED, "");
			return;
		}

		updateStateToClient(JobState.RUNNING, "sorting finished");
	}
}
