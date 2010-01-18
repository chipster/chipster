package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;


public class TsvReadInstructions extends ReadInstructions<Float>{
	public TsvReadInstructions(
			FileDefinition fileDef){
		super(fileDef);
		
		parser = new TsvChunkParser(fileDef);
		chunker = new TsvChunker();
		conciser = new SeqLengthIntensityConciser();
	}
}
