package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;


public class FastaChunkParser extends ChunkParser {
	
	private long rowLength;
	private long titleLength;
	
	File file;

	public FastaChunkParser(FileDefinition fileDef, File file) throws IOException {
		super(fileDef);
		
		this.file = file;
		
		BufferedReader reader = new BufferedReader(new FileReader(file));
		
		titleLength = reader.readLine().length() + 1;
		rowLength = reader.readLine().length() + 1;
	}
	
	public FastaChunkParser(FileDefinition fileDef, long rowLength, long titleLength) {
		super(fileDef);
		
		this.titleLength = titleLength;
		this.rowLength = rowLength;
	}

	@Override
	public float getFloat(long pos) {
		return Float.NaN;
	}

	@Override
	public long getLong(long pos) {
				
		return (chunk.readIndex + pos) * rowLength;
	}

	@Override
	public long getPosition(long readIndex, Content rec) {
			
		return readIndex - chunk.readIndex;
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

		byte[] byteValue = new byte[(int) (rowLength)];
		
		//System.out.println(rowLength + ", " + pos + ", " + chunk.content.length);
		
		System.arraycopy(chunk.content, (int)(pos * rowLength), 
				byteValue, 0, byteValue.length);
		
		return byteValue;
	}
	
	public long getRowLength(){
		return rowLength;
	}

	@Override
	public ChunkParser clone() {
		return new FastaChunkParser(fileDef, rowLength, titleLength);		
	}

	public long getTitleLength() {
		return titleLength;
	}
}
