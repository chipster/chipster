package fi.csc.chipster.tools.ngs;

import java.io.File;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.BedLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.ChromosomeNormaliser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TsvSorter;
import fi.csc.microarray.messaging.JobState;

public class SortBed extends JavaAnalysisJobBase {

	public static final ChromosomeNormaliser CHROMOSOME_NORMALISER = new ChromosomeNormaliser() {

		public String normaliseChromosome(String chromosomeName) {

			// Leave prefix as it is
			
			// Remove postfix, if present
			String SEPARATOR = ".";
			if (chromosomeName.contains(SEPARATOR)) {
				chromosomeName = chromosomeName.substring(0, chromosomeName.indexOf(SEPARATOR));
			}
			
			return chromosomeName;
		}
	};
	
	@Override
	public String getSADL() {
		return 	"TOOL SortBed.java: \"Sort BED\" (Sort a BED file by chromosome and start position.)" + "\n" +
				"INPUT regions.bed: \"BED file\" TYPE GENERIC" + "\n" +
				"OUTPUT sorted.bed: \"Sorted BED file\"" + "\n"; 

	}
	
	
	@Override
	protected void execute() { 
		updateStateToClient(JobState.RUNNING, "sorting");


		try {
			// files
			File inputFile = new File(jobWorkDir, analysis.getInputFiles().get(0).getFileName()); 
			File outputFile = new File(jobWorkDir, analysis.getOutputFiles().get(0).getFileName().getID()); 

			// run sort
			//BEDParser increments coordinates by one, but it's not a problem because only its column order is used
			new TsvSorter().sort(
					inputFile, outputFile, CHROMOSOME_NORMALISER,
					BedLineParser.Column.CHROMOSOME.ordinal(), BedLineParser.Column.START.ordinal(), new BedLineParser(false));

		} catch (Exception e) {
			updateState(JobState.FAILED, e.getMessage());
			return;
		}

		updateStateToClient(JobState.RUNNING, "sorting finished");
	}
}
