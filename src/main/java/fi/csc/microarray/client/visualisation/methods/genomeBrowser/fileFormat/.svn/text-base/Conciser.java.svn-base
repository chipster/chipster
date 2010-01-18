package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.ByteChunk;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;

public interface Conciser<T> {
	public T concise(ByteChunk chunk, ChunkParser parser);
	public Region getRegion(ChunkParser parser, long readIndex);
}
