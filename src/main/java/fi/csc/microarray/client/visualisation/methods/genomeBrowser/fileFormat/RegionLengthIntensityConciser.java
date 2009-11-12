package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.ByteChunk;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;

public class RegionLengthIntensityConciser implements Conciser<Float> {

	public Float concise(ByteChunk chunk, ChunkParser parser) {
			
			Content start = Content.BP_START;
			Content end = Content.BP_END;
			
			long total = 0;
			int i;
			long minBp = -1;
			long maxBp = -1;
			
			parser.setChunk(chunk);
			
			long rowCount = parser.getReadCount();
			
			for (i = 0; i < rowCount; i++){
				
				long startBp = parser.getLong(i + chunk.readIndex, start);
				long length = parser.getLong(i + chunk.readIndex, end) - startBp;

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
		
		return new Region(parser.getLong(readIndex, Content.BP_START),
				parser.getLong(readIndex, Content.BP_END));
	}
}
