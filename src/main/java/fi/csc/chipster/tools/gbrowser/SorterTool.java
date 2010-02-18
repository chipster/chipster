package fi.csc.chipster.tools.gbrowser;

import java.io.File;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ConstantRowLengthParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileDefinition;
import fi.csc.microarray.messaging.JobState;

public class SorterTool extends JavaAnalysisJobBase {

	private ConstantRowLengthParser[] parsers = {
			new ElandParser()
	};
	
	@Override
	public String getVVSADL() {
		
		
		StringBuffer params = new StringBuffer();
		
		for(int i = 0; i < parsers.length; i++) {
			
			params.append(parsers[i].getName());
			
			if (i < parsers.length - 1) {
				params.append(", ");
			}
		}
		
		return 	" ANALYSIS Utils/Sort (Sort primarily using chromosome and secondarily using start " +
				"location of the feature. File format is used to find columns containing " +
				"chromosome and start location. )" + "\n" +
				
				" INPUT GENERIC input.tsv OUTPUT output.tsv" + "\n" +
				" PARAMETER file.format [" + params + "] DEFAULT params[0] ()";
 
	}

	@Override
	protected void execute() throws Exception {
		updateState(JobState.RUNNING, "Sorting file", true);
		
		File inputFile = new File(jobWorkDir, "input.tsv");
		File tmpFile = new File(jobWorkDir, "tmp.tsv");
		File outputFile = new File(jobWorkDir, "output.tsv");		
		
		
		FileDefinition def = null;
		
		for(int i = 0; i < parsers.length; i++) {
			
			if (parsers[i].getName().equals(inputMessage.getParameters().get(0))) {
				def = parsers[i].getFileDefinition();
			}
		}		

		new TsvSorter().sort(inputFile, tmpFile, 
				def.indexOf(ColumnType.CHROMOSOME), def.indexOf(ColumnType.BP_START));
		
		
		
		//TODO Remove following after visualisation is able to read tsv files
		if (inputMessage.getParameters().get(1).equals("yes")) {
			
		
			updateState(JobState.RUNNING, "Converting to constant", true);
			
			//Minus one to skip the last column ćontaining only a new line character
			int[] fieldLengths = new int[def.size() - 1];
			
			for (int i = 0; i < def.size() - 1; i++ ) {
				
				fieldLengths[i] = def.get(i).length;
			}
				
			TsvToConstant.convert(tmpFile, outputFile, fieldLengths);
		}					
		
		updateState(JobState.RUNNING, "Sorting finished", true);		
	}
}
