package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Arrays;


public class ConstantLengthBowtieReadInstructions extends ConstantLengthReadInstructions{
	
	public static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new DataFieldDef[] {

							new DataFieldDef(Content.ID, Type.STRING, 16),
							new DataFieldDef(Content.STRAND, Type.STRING, 2),
							new DataFieldDef(Content.SKIP, Type.STRING, 16),
							new DataFieldDef(Content.BP_START, Type.LONG, 16),
							new DataFieldDef(Content.SEQUENCE, Type.STRING, 50),
							new DataFieldDef(Content.SKIP, Type.STRING, 50),
							new DataFieldDef(Content.SKIP, Type.STRING, 16),
							new DataFieldDef(Content.SKIP, Type.STRING, 16),
							new DataFieldDef(Content.SKIP, Type.NEWLINE, 1)							
					}));

	public ConstantLengthBowtieReadInstructions() {
		super(fileDef);		
	}
}