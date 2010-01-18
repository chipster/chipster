package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;


public class ConstantLengthReadInstructions extends ReadInstructions<Float>{
	public ConstantLengthReadInstructions(FileDefinition fileDef){
		super(fileDef);
		
		parser = new ConstantLengthChunkParser(fileDef);
		chunker = new ConstantLengthChunker((ConstantLengthChunkParser)parser);
		conciser = new SeqLengthIntensityConciser();
	}
}
