package fi.csc.chipster.tools.ngs;

import java.io.File;
import java.io.IOException;

import fi.csc.chipster.tools.gbrowser.ChromosomeNormaliser;
import fi.csc.chipster.tools.gbrowser.SamBamUtils;
import fi.csc.chipster.tools.gbrowser.TsvSorter;
import fi.csc.chipster.tools.gbrowser.SamBamUtils.SamBamUtilState;
import fi.csc.chipster.tools.gbrowser.SamBamUtils.SamBamUtilStateListener;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;

public class LocalNGSPreprocess implements Runnable {

	public static final ChromosomeNormaliser CHROMOSOME_NORMALISER = new ChromosomeNormaliser() {

		public String normaliseChromosome(String chromosomeName) {
			
			// Add prefix, if it is missing
			String CHROMOSOME_NAME_PREFIX = "chr";
			if (!chromosomeName.startsWith(CHROMOSOME_NAME_PREFIX)) {
				chromosomeName = CHROMOSOME_NAME_PREFIX + chromosomeName;
			}
			
			// Remove postfix, if present
			String SEPARATOR = ".";
			if (chromosomeName.contains(SEPARATOR)) {
				chromosomeName = chromosomeName.substring(0, chromosomeName.indexOf(SEPARATOR));
			}
			
			return chromosomeName;
		}
	};

	private static TsvParser[] parsers = {
			new ElandParser()
	};
	
	private Task task;
	
	public LocalNGSPreprocess(Task task) {
		this.task = task;
	}
	
	public static String getSADL() {
		
		StringBuffer fileFormats = new StringBuffer();
		for (int i = 0; i < parsers.length; i++) {
			fileFormats.append(parsers[i].getName() + ": " + parsers[i].getName());
			
			if (i < parsers.length - 1) {
				fileFormats.append(", ");
			}
		}
		
//		String description = "Chipster genome browser is able to show BAM and BED files. BAM files need to be " + 
//		"sorted and indexed, and Chipster can perform this preprocessing for you.\n\n " +
//		"Please note that preprocessing BAM files can take several minutes depending " +
//		"on the file size. Also BED files need to be preprocessed prior to viewing, " +
//		"but this is very quick.";
		String description = "<p>Does my NGS data need to preprocessed by Chipster?</p><br/>" + 

				"<p>-SAM files: yes</p>" +
				"<p>-BAM files: yes, unless your file is already sorted and you have an index file for it</p>" +
				"<p>-BED files: yes, unless your file is already sorted</p>";
		
		
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
					preprocessBed(dataManager, inputFile);
					
				} else if ("bai".equals(extension)) {
					preprocessBai(dataManager, inputFile);
				} else {
					preprocessReads(dataManager, inputFile, extension);
				}
					
				
			}
		} catch (Exception e) {
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
			 
		}, CHROMOSOME_NORMALISER);
		
		if (SamBamUtils.isSamBamExtension(extension)) {
			samBamUtil.preprocessSamBam(inputFile, outputFile, indexOutputFile);
			
		} else {
			// Assume ELAND format
			samBamUtil.preprocessEland(inputFile, outputFile, indexOutputFile);
		}

		// create outputs in the client
		DataBean outputBean = dataManager.createDataBean(outputName, outputFile);
		DataBean indexOutputBean = dataManager.createDataBean(indexOutputName, indexOutputFile);
		
		// create new operation instance, without any inputs FIXME parameters are lost, sucks
		outputBean.setOperation(new Operation(Session.getSession().getApplication().getOperationDefinition(task.getOperationID()), new DataBean[] {}));
		indexOutputBean.setOperation(new Operation(Session.getSession().getApplication().getOperationDefinition(task.getOperationID()), new DataBean[] {}));
		dataManager.getRootFolder().addChild(outputBean);
		dataManager.getRootFolder().addChild(indexOutputBean);
	}
	
	private void preprocessBed(DataManager dataManager, File inputFile) throws Exception {
		
		String outputName = generateFilename(inputFile, "bed");
		File outputFile = dataManager.createNewRepositoryFile(outputName);		

		// Sort
		new TsvSorter().sort(inputFile, outputFile, new BEDParser(), CHROMOSOME_NORMALISER);
		
		// Create outputs in the client
		DataBean outputBean = dataManager.createDataBean(outputName, outputFile);
		
		// Create new operation instance, without any inputs FIXME parameters are lost, sucks
		outputBean.setOperation(new Operation(Session.getSession().getApplication().getOperationDefinition(task.getOperationID()), new DataBean[] {}));
		dataManager.getRootFolder().addChild(outputBean);
	}

	private void preprocessBai(DataManager dataManager, File inputFile) throws Exception {
		String outputName = inputFile.getName();		
		File outputFile = dataManager.createNewRepositoryFile(outputName);		
		
		// Create outputs in the client
		DataBean outputBean = dataManager.createDataBean(outputName, outputFile);
		
		// Create new operation instance, without any inputs FIXME parameters are lost, sucks
		outputBean.setOperation(new Operation(Session.getSession().getApplication().getOperationDefinition(task.getOperationID()), new DataBean[] {}));
		dataManager.getRootFolder().addChild(outputBean);
	}
	
}
