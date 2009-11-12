package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Arrays;


public class ConstantRefGeneInstructions extends ConstantLengthReadInstructions{
	
	public ConstantRefGeneInstructions(){
		super(fileDef);

		conciser = new RegionLengthIntensityConciser();
	}

	public static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new DataFieldDef[] {

							new DataFieldDef(Content.CHROMOSOME, Type.STRING, 16),
							new DataFieldDef(Content.SKIP, Type.STRING, 16),
							new DataFieldDef(Content.DESCRIPTION, Type.STRING, 16),
							new DataFieldDef(Content.BP_START, Type.LONG, 16),
							new DataFieldDef(Content.BP_END, Type.LONG, 16),
							new DataFieldDef(Content.SKIP, Type.FLOAT, 16),
							new DataFieldDef(Content.STRAND, Type.STRING, 2),
							new DataFieldDef(Content.SKIP, Type.STRING, 2),
							new DataFieldDef(Content.ID, Type.STRING, 64),
							new DataFieldDef(Content.SKIP, Type.NEWLINE, 1)
					}));
}