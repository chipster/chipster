package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.Iterator;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher.ByteChunk;

public class TsvChunkParser extends ChunkParser {

	public static final long NEXT = -1;
	private String[] line;
	private Iterator<String> lineIterator;
	private BufferedReader input;

	public TsvChunkParser(FileDefinition fileDef) {
		super(fileDef);
	}

	@Override
	public float getFloat(long pos) {
		return new Float(getString(pos));
	}

	@Override
	public long getLong(long pos) {
		return new Long(getString(pos));
	}

	@Override
	public long getPosition(long readIndex, Content rec) {
		if (readIndex == NEXT) {
			return fileDef.indexOf(rec);
		} else {
			throw new UnsupportedOperationException(
					"Only iterative use supported");
		}
	}

	public long getNextFieldPosition() {
		return NEXT;
	}

	@Override
	public long getReadCount() {
		throw new UnsupportedOperationException("Only iterative use supported");
	}

	@Override
	public String getString(long pos) {
		if (pos == NEXT) {
			return getNextString();
		} else {
			return getField(pos);
		}
	}
	
	private String getField(long pos){
		if(line == null){
			readNextLine();
		}
		
		return line[(int) pos];
	}

	private String getNextString() {

		if (lineIterator == null || !lineIterator.hasNext()) {
			readNextLine();
		}

		return lineIterator.next();
	}

	public boolean readNextLine() {
		if (input == null) {
			input = new BufferedReader(new InputStreamReader(
					new ByteArrayInputStream(chunk.content)));
		}
		String line;
		try {

			if ((line = input.readLine()) != null) {
				
				this.line = line.split("\t");
				lineIterator = Arrays.asList(this.line).iterator();
				return true;
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
		return false;
	}
	
	@Override
	public void setChunk(ByteChunk chunk){
		super.setChunk(chunk);
		
		this.input = null;
		this.line = null;
		this.lineIterator = null;
	}

	@Override
	public ChunkParser clone() {
		return new TsvChunkParser(fileDef);
	}
}
