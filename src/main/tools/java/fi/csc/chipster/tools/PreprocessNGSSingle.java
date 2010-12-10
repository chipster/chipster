package fi.csc.chipster.tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import fi.csc.chipster.tools.gbrowser.SamBamUtils;
import fi.csc.chipster.tools.gbrowser.TsvSorter;
import fi.csc.chipster.tools.gbrowser.SamBamUtils.SamBamUtilState;
import fi.csc.chipster.tools.gbrowser.SamBamUtils.SamBamUtilStateListener;
import fi.csc.microarray.analyser.AnalysisDescription.InputDescription;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.IOUtils;

public class PreprocessNGSSingle extends JavaAnalysisJobBase {

	@Override
	public String getSADL() {
		return 	"TOOL \"Preprocess\" / PreprocessNGSSingle.java: \"Preprocess NGS\" (Preprocess NGS data to be visualised in the genome browser.)" + "\n" +
				"INPUT data.bam: \"Input bam file\" TYPE GENERIC" + "\n" +
				"OUTPUT preprocessed.bam: \"Preprocessed bam file\"" + "\n" + 
		        "OUTPUT preprocessed.bam.bai: \"Preprocessed bam index file\"" + "\n"; 

	}

	
	
	@Override
	protected void execute() { 
		updateState(JobState.RUNNING, "Sorting file");
		
//
//		for (InputDescription input: analysis.getInputFiles()) {
//			File inputFile = new File(jobWorkDir, input.getFileName()); 
//			File outputFile = new File(jobWorkDir, input.getFileName().substring("in-".length()));		
//
//			
//			// run sorter
//			try {
//				new TsvSorter().sort(inputFile, outputFile, parser);
//
//			
//			
//				String indexOutputName = outputName + ".bai";
//
//				File indexOutputFile = dataManager.createNewRepositoryFile(indexOutputName);
//
//				// Run preprocessing
//				SamBamUtils samBamUtil= new SamBamUtils(new SamBamUtilStateListener() {
//
//					@Override
//					public void stateChanged(SamBamUtilState newState) {
//						task.setStateDetail(newState.getState() + " " + newState.getPercentage());
//					}
//					 
//				}, CHROMOSOME_NORMALISER);
//				
//				if (SamBamUtils.isSamBamExtension(extension)) {
//					samBamUtil.preprocessSamBam(inputFile, outputFile, indexOutputFile);
//					
//				} else {
//					// Assume ELAND format
//					samBamUtil.preprocessEland(inputFile, outputFile, indexOutputFile);
//				}
//
//				// create outputs in the client
//				DataBean outputBean = dataManager.createDataBean(outputName, outputFile);
//				DataBean indexOutputBean = dataManager.createDataBean(indexOutputName, indexOutputFile);
//				
//				// create new operation instance, without any inputs FIXME parameters are lost, sucks
//				outputBean.setOperation(new Operation(Session.getSession().getApplication().getOperationDefinition(task.getOperationID()), new DataBean[] {}));
//				indexOutputBean.setOperation(new Operation(Session.getSession().getApplication().getOperationDefinition(task.getOperationID()), new DataBean[] {}));
//				dataManager.getRootFolder().addChild(outputBean);
//				dataManager.getRootFolder().addChild(indexOutputBean);
//			
//			
//			
//			
//			
//			
//			
//			
			
			
			
//			
//			
//			
//			} catch (Exception e) {
//				updateState(JobState.FAILED, e.getMessage());
//				return;
//			}
//		}
//		
//		
//		
//		updateState(JobState.RUNNING, "sort finished");
	}
}
