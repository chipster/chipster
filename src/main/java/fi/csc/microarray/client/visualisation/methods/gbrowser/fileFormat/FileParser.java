package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Collection;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ByteChunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RowRegion;

public abstract class FileParser {

	protected ByteChunk chunk;

	public abstract long getRowIndex(long filePosition);
	public abstract long getFilePosition(long readIndex);
	public abstract long getChunkRowCount();
	public abstract Object get(long readIndex, ColumnType field);
	public abstract Map<ColumnType, Object> getValues(long readIndex, Collection<ColumnType> requestedFields);
	public abstract RegionContent[] concise(BpCoordRegion readIndexRegion);
	public abstract int getChunkMaxByteLength();
	public abstract BpCoordRegion getBpRegion(long readIndex);
	public abstract int getRowByteLength();
	public abstract String getName();
	public abstract FileParser clone();

	public void setChunk(ByteChunk chunk) {
		this.chunk = chunk;
	}

	public long getRowIndex() {
		return chunk.rowIndex;
	}

	public RowRegion getChunkRegionMiddleOf(RowRegion rowIndexes) {

		RowRegion nodeRows = new RowRegion();

		// round to next chunk split
		nodeRows.start = (long) Math.ceil(rowIndexes.getMid() / getChunkMaxRowCount()) * getChunkMaxRowCount();
		nodeRows.end = nodeRows.start + getChunkMaxRowCount() - 1;

		return nodeRows;
	}

	public long getChunkMaxRowCount() {
		return getChunkMaxByteLength() / getRowByteLength();
	}


	public int getChildCount(RowRegion subtreeReadIndexes) {

		if (subtreeReadIndexes.getLength() <= getChunkMaxByteLength() / getRowByteLength()) {
			return 0;
		} else if (subtreeReadIndexes.getLength() <= getChunkMaxRowCount() * 2) {
			return 1;
		} else {
			return 2;
		}
	}

}
