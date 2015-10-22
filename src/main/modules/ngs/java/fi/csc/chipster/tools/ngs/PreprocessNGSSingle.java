package fi.csc.chipster.tools.ngs;

import java.io.File;

import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils.SamBamUtilState;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils.SamBamUtilStateListener;
import fi.csc.microarray.comp.java.JavaCompJobBase;
import fi.csc.microarray.messaging.JobState;

public class PreprocessNGSSingle extends JavaCompJobBase {

	@Override
	public String getSADL() {
		return 	"TOOL PreprocessNGSSingle.java: \"Convert SAM to BAM, sort and index BAM\" (Converts SAM to BAM and sorts and indexes BAM files. Please note that this preprocessing is required for visualising the data in the Chipster Genome browser. This tools is based on the Picard package.)" + "\n" +
				"INPUT data.bam: \"Input bam file\" TYPE GENERIC" + "\n" +
				"OUTPUT preprocessed.bam: \"Preprocessed bam file\"" + "\n" + 
		        "OUTPUT preprocessed.bam.bai: \"Preprocessed bam index file\"" + "\n"; 
	}
	
	
	@Override
	protected void execute() { 
		updateStateToClient(JobState.RUNNING, "preprocessing");


		try {
			// files
			File inputFile = new File(jobWorkDir, toolDescription.getInputFiles().get(0).getFileName()); 
			File outputFile = new File(jobWorkDir, toolDescription.getOutputFiles().get(0).getFileName().getID()); 
			File indexOutputFile = new File(jobWorkDir, toolDescription.getOutputFiles().get(1).getFileName().getID());


			// run preprocessing
			SamBamUtils samBamUtil= new SamBamUtils(new SamBamUtilStateListener() {

				@Override
				public void stateChanged(SamBamUtilState newState) {
					// update detail state
					updateStateDetailToClient("preprocess: " + newState.getState());
				}

			});

			samBamUtil.preprocessSamBam(inputFile, outputFile, indexOutputFile);


		} catch (Exception e) {
			updateState(JobState.FAILED, e.getMessage());
			return;
		}

		updateStateToClient(JobState.RUNNING, "preprocessing finished");
	}
}
