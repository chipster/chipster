package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Collection;
import java.util.Map;

public abstract class TsvFileParser extends FileParser{
	
	public FileDefinition fileDef;
	
	@Override
	public Object get(long readIndex, ColumnType field) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public long getChunkRowCount() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getRowByteLength() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public Map<ColumnType, Object> getValues(long readIndex,
			Collection<ColumnType> requestedFields) {
		// TODO Auto-generated method stub
		return null;
	}
}
