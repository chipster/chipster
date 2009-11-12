package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Arrays;

public class FastaFsfReadInstructions extends ConstantLengthReadInstructions{
	
	//int[] fieldLengths = new int[] { 32, 64, 8, 8, 8, 8, 16, 16, 2};
	
	public static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new DataFieldDef[] {

							new DataFieldDef(Content.SEQUENCE, Type.STRING, 64),
							new DataFieldDef(Content.BP_START, Type.LONG, 16),							
							new DataFieldDef(Content.SKIP, Type.NEWLINE, 1),
														
			}));

	public FastaFsfReadInstructions() {
		super(fileDef);		
	}
}