package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Collection;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.ByteChunk;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RowRegion;

public abstract class FileParser {
	
	protected ByteChunk chunk;	
	
	public abstract long getRowIndex(long filePosition);
	
	public abstract long getFilePosition(long readIndex);
		
	public abstract long getChunkRowCount();
	
	public abstract Object get(long readIndex, ColumnType field);

	public abstract Map<ColumnType, Object> getValues(long readIndex,
			Collection<ColumnType> requestedFields);
	
	public abstract RegionContent[] concise(BpCoordRegion readIndexRegion);

	public abstract int getChunkMaxByteLength();

	public abstract BpCoordRegion getBpRegion(long readIndex);
	
	public void setChunk(ByteChunk chunk) {
		this.chunk = chunk;
	}
	
	public abstract FileParser clone();
	
	public long getRowIndex(){
		return chunk.rowIndex;
	}
	
	public RowRegion getChunkRegionMiddleOf(RowRegion rowIndexes){

		RowRegion nodeRows = new RowRegion();
		
		//Round to next chunk split
		nodeRows.start = 
			(long)Math.ceil(rowIndexes.getMid() / getChunkMaxRowCount()) * getChunkMaxRowCount();
		nodeRows.end = nodeRows.start + getChunkMaxRowCount() - 1;

		return nodeRows;
	}
	
	public long getChunkMaxRowCount() { 
		return getChunkMaxByteLength() / getRowByteLength();
	}

	public abstract int getRowByteLength();

	public int getChildCount(RowRegion subtreeReadIndexes){

		if(subtreeReadIndexes.getLength() <= getChunkMaxByteLength() / getRowByteLength()){
			return 0;
		} else if(subtreeReadIndexes.getLength() <= getChunkMaxRowCount() * 2){
			return 1;
		} else {
			return 2;
		}
	}
	
	public abstract String getName();
}
