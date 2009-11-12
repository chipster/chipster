package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;


public class FastaReadInstructions extends TsvReadInstructions {

	public static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new DataFieldDef[] {

							new DataFieldDef(Content.SEQUENCE, Type.STRING)						
			}));

	

	public FastaReadInstructions(File file) throws IOException {
		super(fileDef);		
		
		parser = new FastaChunkParser(fileDef, file);
		chunker = new FastaChunker((FastaChunkParser) parser);
		conciser = new DummyConciser();
	}
}

