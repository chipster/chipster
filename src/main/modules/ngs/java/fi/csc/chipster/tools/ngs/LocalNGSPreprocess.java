package fi.csc.chipster.tools.ngs;

import java.io.File;
import java.io.IOException;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.AbstractTsvLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.BedLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.GtfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.VcfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils.SamBamUtilState;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils.SamBamUtilStateListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.TsvSorter;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;

public class LocalNGSPreprocess implements Runnable {
	
	private static final Logger logger = Logger.getLogger(LocalNGSPreprocess.class);
	
	private Task task;
	
	public LocalNGSPreprocess(Task task) {
		this.task = task;
	}
	
	public static String getSADL() {
		
//		String description = "Chipster genome browser is able to show BAM and BED files. BAM files need to be " + 
//		"sorted and indexed, and Chipster can perform this preprocessing for you.\n\n " +
//		"Please note that preprocessing BAM files can take several minutes depending " +
//		"on the file size. Also BED files need to be preprocessed prior to viewing, " +
//		"but this is very quick.";
		String description = "<p>Does my NGS data need to preprocessed by Chipster?</p><br/>" + 

				"<p>-SAM files: yes</p>" +
				"<p>-BAM files: yes, unless your file is already sorted and you have an index file for it</p>" +
				"<p>-BED, GTF and VCF files: yes, unless your file is already sorted</p>";
		
		
		return 	"TOOL LocalNGSPreprocess.java: \"NGS Preprocess\" (" + description + ")" + "\n" +
				"INPUT input{...}.txt: \"Input NGS data\" TYPE GENERIC" + "\n" +
				"OUTPUT ngs-preprocess.txt: \"Preprocessed NGS data\"" + "\n" +
				"OUTPUT phenodata.tsv: \"Phenodata\"";
 	}

	public void run() {
		DataManager dataManager = Session.getSession().getDataManager();
		
		try {

			for (DataBean inputDataBean : task.getInputs()) {
				File inputFile = dataManager.getLocalFile(inputDataBean);
				String extension = inputFile.getName().substring(inputFile.getName().lastIndexOf(".") + 1);
				
				if ("bed".equals(extension)) {
					preprocess(dataManager, inputFile, "bed", new BedLineParser(false), 
							BedLineParser.Column.CHROM.ordinal(), BedLineParser.Column.CHROM_START.ordinal());
					
				} else if ("gtf".equals(extension)) {
					preprocess(dataManager, inputFile, "gtf", new GtfLineParser(), 
							GtfLineParser.Column.SEQNAME.ordinal(), GtfLineParser.Column.START.ordinal());
					
				} else if ("vcf".equals(extension)) {
					preprocess(dataManager, inputFile, "vcf", new VcfLineParser(), 
							VcfLineParser.Column.CHROM.ordinal(), VcfLineParser.Column.POS.ordinal());
					
				} else if ("bai".equals(extension)) {
					preprocessBai(dataManager, inputFile);
				} else {
					preprocessReads(dataManager, inputFile, extension);
				}
					
				
			}
		} catch (Exception e) {
			logger.error(e);
			throw new RuntimeException(e);
		}
	}

	private String generateFilename(File inputFile, String extension) {
		
		// Strip extension and add "-preprocessed"
		String outputName;
		int fileExtensionStartPosition = inputFile.getName().lastIndexOf(".");
		if (fileExtensionStartPosition != -1) {
			outputName = inputFile.getName().substring(0, fileExtensionStartPosition) + "-preprocessed";
		} else {
			outputName = inputFile.getName() + "-preprocessed";
		}
		
		// Replace new extension
		outputName = outputName + "." + extension;

		return outputName;
	}

	private void preprocessReads(DataManager dataManager, File inputFile, String extension) throws IOException, MicroarrayException {
		
		String outputName = generateFilename(inputFile, "bam");
		String indexOutputName = outputName + ".bai";

		File outputFile = dataManager.createNewRepositoryFile(outputName);		
		File indexOutputFile = dataManager.createNewRepositoryFile(indexOutputName);

		// Run preprocessing
		SamBamUtils samBamUtil= new SamBamUtils(new SamBamUtilStateListener() {

			@Override
			public void stateChanged(SamBamUtilState newState) {
				task.setStateDetail(newState.getState() + " " + newState.getPercentage());
			}
			 
		});
		
		if (SamBamUtils.isSamBamExtension(extension)) {
			samBamUtil.preprocessSamBam(inputFile, outputFile, indexOutputFile);
			
		} else {
			// Assume ELAND format
			samBamUtil.preprocessEland(inputFile, outputFile, indexOutputFile);
		}

		// create outputs in the client
		DataBean outputBean = dataManager.createDataBean(outputName, outputFile);
		DataBean indexOutputBean = dataManager.createDataBean(indexOutputName, indexOutputFile);
		
		// Create new operation instance, without any inputs FIXME parameters are lost, sucks create OperationRecord directly
		OperationRecord operationRecord = new OperationRecord(new Operation(Session.getSession().getApplication().getOperationDefinition(task.getOperationID()), new DataBean[] {}));
		outputBean.setOperationRecord(operationRecord);
		indexOutputBean.setOperationRecord(operationRecord);
		//Chipster2 backport fix
		dataManager.getRootFolder().addChild(outputBean);
		dataManager.getRootFolder().addChild(indexOutputBean);
	}
	
	private void preprocess(DataManager dataManager, File inputFile, String fileExtension, AbstractTsvLineParser lineParser, int chrColumn, int startColumn) throws Exception {
		
		String outputName = generateFilename(inputFile, fileExtension);
		File outputFile = dataManager.createNewRepositoryFile(outputName);		

		// Sort
		new TsvSorter().sort(
				inputFile, outputFile, chrColumn, startColumn, lineParser);
		
		createOutput(dataManager, outputName, outputFile);
	}

	private void preprocessBai(DataManager dataManager, File inputFile) throws Exception {
		String outputName = inputFile.getName();		
		File outputFile = dataManager.createNewRepositoryFile(outputName);		
		
		createOutput(dataManager, outputName, outputFile);
	}
	
	private void createOutput(DataManager dataManager, String outputName,
			File outputFile) throws MicroarrayException {
		// Create outputs in the client
		DataBean outputBean = dataManager.createDataBean(outputName, outputFile);
		
		// Create new operation instance, without any inputs FIXME parameters are lost, sucks create OperationRecord directly
		outputBean.setOperationRecord(new OperationRecord(new Operation(Session.getSession().getApplication().getOperationDefinition(task.getOperationID()), new DataBean[] {})));
		//Chipster2 backport fix
		dataManager.getRootFolder().addChild(outputBean);
	}
}
