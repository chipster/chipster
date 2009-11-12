package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Arrays;

public class ElandReadInstructions extends ConstantLengthReadInstructions{
	
	//int[] fieldLengths = new int[] { 32, 64, 8, 8, 8, 8, 16, 16, 2};
	
	public static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new DataFieldDef[] {

							new DataFieldDef(Content.ID, Type.STRING, 32),
							new DataFieldDef(Content.SEQUENCE, Type.STRING, 64),
							new DataFieldDef(Content.QUALITY, Type.STRING, 8),
							new DataFieldDef(Content.SKIP, Type.STRING, 8),
							new DataFieldDef(Content.SKIP, Type.STRING, 8),
							new DataFieldDef(Content.SKIP, Type.STRING, 8),
							new DataFieldDef(Content.CHROMOSOME, Type.STRING, 16),
							new DataFieldDef(Content.BP_START, Type.LONG, 16),							
							new DataFieldDef(Content.STRAND, Type.STRING, 2),
							new DataFieldDef(Content.SKIP, Type.NEWLINE, 1),
														
			}));

	public ElandReadInstructions() {
		super(fileDef);		
	}
}