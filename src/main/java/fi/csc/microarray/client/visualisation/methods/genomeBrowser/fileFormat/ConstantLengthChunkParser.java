package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;


public class ConstantLengthChunkParser extends ChunkParser {
	
	private long rowLength;

	public ConstantLengthChunkParser(FileDefinition fileDef) {
		super(fileDef);
		
		for(DataFieldDef col : fileDef){
			col.offset = rowLength;
			rowLength += col.length;
		}
	}

	@Override
	public float getFloat(long pos) {
		return new Float(getString(pos));
	}

	@Override
	public long getLong(long pos) {
		
		String str = getString(pos).trim();
		if(!str.equals(""))
			return new Long(str);
		else {
			return Long.MIN_VALUE;
		}
	}

	@Override
	public long getPosition(long readIndex, Content rec) {
		return readIndex * fileDef.size() + fileDef.indexOf(rec);
	}

	@Override
	public long getReadCount() {
		return chunk.length / rowLength; 
	}

	@Override
	public String getString(long pos) {
		return new String(getBytes(pos));
	}
	
	private byte[] getBytes(long pos){
		long row = pos / fileDef.size() - chunk.readIndex;
		int colIndex = (int)(pos % fileDef.size());
		
		byte[] byteValue = new byte[fileDef.get(colIndex).length];
		
		System.arraycopy(chunk.content, (int)(row*rowLength + fileDef.get(colIndex).offset), 
				byteValue, 0, byteValue.length);
		return byteValue;
	}
	
	public long getRowLength(){
		return rowLength;
	}

	@Override
	public ChunkParser clone() {
		return new ConstantLengthChunkParser(fileDef);
	}
}
