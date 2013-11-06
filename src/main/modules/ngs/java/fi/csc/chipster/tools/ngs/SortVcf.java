package fi.csc.chipster.tools.ngs;

import java.io.File;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.VcfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TsvSorter;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.Exceptions;

public class SortVcf extends JavaAnalysisJobBase {
	
	@Override
	public String getSADL() {
		return 	"TOOL SortVcf.java: \"Sort VCF\" (Sort a VCF file by chromosome and position.)" + "\n" +
				"INPUT unsorted.vcf: \"VCF file\" TYPE GENERIC" + "\n" +
				"OUTPUT sorted.vcf: \"Sorted VCF file\"" + "\n"; 

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
					VcfLineParser.Column.CHROM.ordinal(), 
					VcfLineParser.Column.POS.ordinal(), new VcfLineParser());

		} catch (Exception e) {
			getResultMessage().setErrorMessage(Exceptions.getStackTrace(e));
			updateState(JobState.FAILED, "");
			return;
		}

		updateStateToClient(JobState.RUNNING, "sorting finished");
	}
}

