package fi.csc.chipster.tools;

import java.io.File;

import fi.csc.chipster.tools.gbrowser.TsvSorter;
import fi.csc.chipster.tools.gbrowser.TsvToConstant;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ConstantRowLengthParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileDefinition;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.JobState;

public class SorterTool extends JavaAnalysisJobBase {

	private ConstantRowLengthParser[] parsers = {
			new ElandParser()
	};
	
	@Override
	public String getVVSADL() {
		
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
				" PARAMETER file.format [" + fileFormats + "] DEFAULT " + parsers[0].getName() + " ()" + "\n" +
				" PARAMETER add.spaces [yes, no] DEFAULT yes ()";
 	}

	@Override
	protected void execute() { 
		updateState(JobState.RUNNING, "Sorting file", true);
		
		File inputFile = new File(jobWorkDir, "input.tsv");
		File tmpFile = new File(jobWorkDir, "tmp.tsv");
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
			new TsvSorter().sort(inputFile, tmpFile, 
					def.indexOf(ColumnType.CHROMOSOME), def.indexOf(ColumnType.BP_START));
		} catch (Exception e) {
			updateState(JobState.FAILED, e.getMessage(), true);
		}
		
		// convert tabs to constant length spaces 
		// TODO remove after the visualization is able to read tsv files
		if (inputMessage.getParameters().get(1).equals("yes")) {
		
			updateState(JobState.RUNNING, "converting tabs to constant white space", true);
			
			// get the field lengths
			// minus one to skip the last column ćontaining only a new line character
			int[] fieldLengths = new int[def.size() - 1];
			for (int i = 0; i < def.size() - 1; i++ ) {
				fieldLengths[i] = def.get(i).length;
			}
				
			// convert
			try {
				TsvToConstant.convert(tmpFile, outputFile, fieldLengths);
			} catch (MicroarrayException e) {
				updateState(JobState.FAILED, e.getMessage(), true);
			}
		} 
		
		// don't convert tabs, try to move the tmp file to the output file
		else {
			if (tmpFile.renameTo(outputFile)) {
				updateState(JobState.RUNNING, "sort finished", true);
			} else {
				updateState(JobState.FAILED, "sort failed", true);
			}
		}
	}
}
