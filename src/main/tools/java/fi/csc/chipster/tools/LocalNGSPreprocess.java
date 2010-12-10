package fi.csc.chipster.tools;

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
		
		return 	"TOOL \"Preprocess\" / LocalNGSPreprocess.java: \"NGS Preprocess\" (Sort primarily using chromosome and secondarily using start " +
				"location of the feature. File format is used to find columns containing " +
				"chromosome and start location. )" + "\n" +
				"INPUT input{...}.txt: \"Input NGS data\" TYPE GENERIC" + "\n" +
				"OUTPUT ngs-preprocess.txt: \"Preprocessed NGS data\"" + "\n" +
				"OUTPUT phenodata.tsv: \"Phenodata\"";
 	}

	public void run() {
		DataManager dataManager = Session.getSession().getDataManager();
		
		try {

			for (DataBean inputDataBean : task.getInputs()) {
				File inputFile = dataManager.getLocalFile(inputDataBean);
				String extension = inputFile.getName().substring(inputFile.getName().lastIndexOf("."));
				
				if (".bed".equals(extension)) {
					preprocessBed(dataManager, inputFile);
					
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

}
