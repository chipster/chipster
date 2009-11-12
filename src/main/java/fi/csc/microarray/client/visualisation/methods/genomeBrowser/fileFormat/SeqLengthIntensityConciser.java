package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.ByteChunk;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;

public class SeqLengthIntensityConciser implements Conciser<Float> {

	public Float concise(ByteChunk chunk, ChunkParser parser) {
			
			Content start = Content.BP_START;
			Content seq = Content.SEQUENCE;
			
			long total = 0;
			int i;
			long minBp = -1;
			long maxBp = -1;
			
			parser.setChunk(chunk);
			
			long rowCount = parser.getReadCount();
			
			for (i = 0; i < rowCount; i++){
				
				long startBp = parser.getLong(i + chunk.readIndex, start);
				long length = parser.getString(i + chunk.readIndex, seq).trim().length();

				total += length;
				if(i == 0){
					minBp = startBp;
				}
				if(i == rowCount - 1){
					maxBp = startBp;			
				}
			}
			
			return total / (float)(maxBp - minBp);
		}
	
	public Region getRegion(ChunkParser parser, long readIndex) {
		
		long startBp = parser.getLong(readIndex, Content.BP_START);
		long length = parser.getString(readIndex, Content.SEQUENCE).trim().length();
		
		return new Region(startBp, startBp + length);
	}
}
