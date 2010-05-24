package fi.csc.chipster.tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import fi.csc.chipster.tools.gbrowser.TsvSorter;
import fi.csc.microarray.analyser.AnalysisDescription.InputDescription;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileDefinition;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.util.IOUtils;

public class PreprocessNGSSingle extends JavaAnalysisJobBase {

	private TsvParser[] parsers = {
			new ElandParser()
	};
	
	@Override
	public String getSADL() {
		
		StringBuffer fileFormats = new StringBuffer();
		for (int i = 0; i < parsers.length; i++) {
			fileFormats.append(parsers[i].getName() + ": " + parsers[i].getName());
			
			if (i < parsers.length - 1) {
				fileFormats.append(", ");
			}
		}
		
		// TODO more verbose name, name of the second parameter
		return 	"TOOL \"Preprocess\" / PreprocessNGSSingle.java: \"Preprocess NGS, treatment only\" (Sort primarily using chromosome and secondarily using start " +
				"location of the feature. File format is used to find columns containing " +
				"chromosome and start location. )" + "\n" +
				"INPUT in-treatment.tsv: \"Treatment\" TYPE GENERIC" + "\n" +
				"INPUT in-control.tsv: \"Control\" TYPE GENERIC" + "\n" +
				"OUTPUT treatment.tsv: \"Treatment\"" + "\n" +
				"OUTPUT control.tsv: \"Control\"" + "\n" +
				"OUTPUT phenodata.tsv: \"Phenodata\"" + "\n" +
				"PARAMETER file.format: \"Data format\" TYPE [" + fileFormats + "] DEFAULT " + parsers[0].getName() + " (Format of the data)" + "\n";
 	}

	
	
	@Override
	protected void execute() { 
		updateState(JobState.RUNNING, "Sorting file");
		
		// get the file format and definitions
		FileDefinition def = null;
		for (int i = 0; i < parsers.length; i++) {
			if (parsers[i].getName().equals(inputMessage.getParameters().get(0))) {
				def = parsers[i].getFileDefinition();
			}
		}		

		// sort
		// FIXME check for optionality
		for (InputDescription input: analysis.getInputFiles()) {
			File inputFile = new File(jobWorkDir, input.getFileName()); 
			File outputFile = new File(jobWorkDir, input.getFileName().substring("in-".length()));		
		
			// run sorter
			try {
				new TsvSorter().sort(inputFile, outputFile, 
						def.indexOf(ColumnType.CHROMOSOME), def.indexOf(ColumnType.BP_START));
			} catch (Exception e) {
				updateState(JobState.FAILED, e.getMessage());
				return;
			}
		}
		
		// generate phenodata
		File phenodataFile = new File(jobWorkDir, "phenodata.tsv");
		FileWriter writer = null;
		try {
			writer = new FileWriter(phenodataFile);
			writer.write("sample" + "\t" + "chiptype" + "\t" + "group" + "\t" + "description" + "\n");
//			for (InputDescription input: analysis.getInputFiles()) {
//				writer.write(input.getFileName() + "\t" + "\n");
//			}
			writer.write("treatment.tsv" + "\t\t" + "treatment" + "\n");
			writer.write("control.tsv" + "\t\t" + "control" + "\n");

			
		} catch (IOException e) {
			updateState(JobState.FAILED, e.getMessage());
			return;
		} finally {
			IOUtils.closeIfPossible(writer);
		}
		
		
		updateState(JobState.RUNNING, "sort finished");
	}
}
