package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Arrays;


public class CytobandReadInstructions extends TsvReadInstructions{

	public static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new DataFieldDef[] {

							new DataFieldDef(Content.CHROMOSOME, Type.STRING),
							new DataFieldDef(Content.BP_START, Type.LONG),
							new DataFieldDef(Content.BP_END, Type.LONG),
							new DataFieldDef(Content.ID, Type.STRING),
							new DataFieldDef(Content.VALUE, Type.STRING)
					}));

	public CytobandReadInstructions() {
		super(fileDef);		
	}
}