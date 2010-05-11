package fi.csc.chipster.tools;

import java.io.File;

import fi.csc.chipster.tools.gbrowser.TsvSorter;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileDefinition;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.messaging.JobState;

public class SorterTool extends JavaAnalysisJobBase {

	private TsvParser[] parsers = {
			new ElandParser()
	};
	
	@Override
	public String getSADL() {
		
		StringBuffer fileFormats = new StringBuffer();
		for (int i = 0; i < parsers.length; i++) {
			fileFormats.append(parsers[i].getName());
			
			if (i < parsers.length - 1) {
				fileFormats.append(", ");
			}
		}
		
		// TODO more verbose name, name of the second parameter
		return 	" ANALYSIS Utils/Sort (Sort primarily using chromosome and secondarily using start " +
				"location of the feature. File format is used to find columns containing " +
				"chromosome and start location. )" + "\n" +
				
				" INPUT GENERIC input.tsv OUTPUT output.tsv" + "\n" +
				" PARAMETER file.format [" + fileFormats + "] DEFAULT " + parsers[0].getName() + " ()" + "\n";
 	}

	@Override
	protected void execute() { 
		updateState(JobState.RUNNING, "Sorting file");
		
		File inputFile = new File(jobWorkDir, "input.tsv");
		File outputFile = new File(jobWorkDir, "output.tsv");		
		
		// get the file format and definitions
		FileDefinition def = null;
		for (int i = 0; i < parsers.length; i++) {
			if (parsers[i].getName().equals(inputMessage.getParameters().get(0))) {
				def = parsers[i].getFileDefinition();
			}
		}		

		// run sorter
		try {
			new TsvSorter().sort(inputFile, outputFile, 
					def.indexOf(ColumnType.CHROMOSOME), def.indexOf(ColumnType.BP_START));
		} catch (Exception e) {
			updateState(JobState.FAILED, e.getMessage());
		}
		
		updateState(JobState.RUNNING, "sort finished");
	}
}
