package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Arrays;

public class BedInstructions extends ConstantLengthReadInstructions{

	public static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new DataFieldDef[] {

							new DataFieldDef(Content.CHROMOSOME, Type.STRING, 16),
							new DataFieldDef(Content.BP_START, Type.LONG, 16),
							new DataFieldDef(Content.BP_END, Type.LONG, 16),
							new DataFieldDef(Content.SKIP, Type.STRING, 64),
							new DataFieldDef(Content.SKIP, Type.STRING, 16),
							new DataFieldDef(Content.SKIP, Type.STRING, 16),
							new DataFieldDef(Content.SKIP, Type.NEWLINE, 1)
							
			}));
	
	public BedInstructions() {
		super(fileDef);		
		conciser = new RegionLengthIntensityConciser();
	}
}


