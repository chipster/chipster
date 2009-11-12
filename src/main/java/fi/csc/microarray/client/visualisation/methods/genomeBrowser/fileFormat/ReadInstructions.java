package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;


public class ReadInstructions<T> {
	protected ChunkParser parser;
	public Chunker chunker;
	public Conciser<T> conciser;
	public FileDefinition fileDef;	
	
	public ReadInstructions(FileDefinition fileDef){
		
		this.fileDef = fileDef;
	}
	
	
	/**
	 * Other objects are only read or have no attributes and are thus thread safe
	 */
	public ChunkParser getParser(){
		return (ChunkParser)(parser.clone());
	}
}
