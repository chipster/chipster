package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public class FastaParser extends ConstantRowLengthParser{

	private static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new ColumnDefinition[] {

							new ColumnDefinition(ColumnType.SEQUENCE, Type.STRING, 64),
							new ColumnDefinition(ColumnType.BP_START, Type.LONG, 16),							
							new ColumnDefinition(ColumnType.SKIP, Type.NEWLINE, 1),

					}));

	private Integer rowLength;
	private long titleLength;
	private Chromosome chr;
	
	public FastaParser(File file) throws IOException {
			super(fileDef);
			
			BufferedReader reader = new BufferedReader(new FileReader(file));
			
			titleLength = reader.readLine().length() + 1;
			rowLength = reader.readLine().length() + 1;
	}
	
	public FastaParser(long titleLength, Integer rowLength, Chromosome chr) {
		super(fileDef);
		
		this.titleLength = titleLength;
		this.rowLength = rowLength;
		this.chr = chr;
}

	@Override
	public int getChunkMaxByteLength() {
		return (int)getRowByteLength() * 128;
	}

	@Override
	public long getFilePosition(long readIndex) {
		return readIndex * rowLength + titleLength;
	}

	@Override
	public long getRowIndex(long filePosition) {
		return (filePosition - titleLength) / rowLength;
	}

	@Override
	public RegionContent[] concise(BpCoordRegion readIndexRegion) {

		//Return empty table, otherwise TreeNode gets stuck in calling this again
		return new RegionContent[0];
	}

	@Override
	public BpCoordRegion getBpRegion(long rowIndex) {

		long startBp = (Long)get(rowIndex, ColumnType.BP_START);
		long length = ((String)get(rowIndex, ColumnType.SEQUENCE)).trim().length();

		return new BpCoordRegion(startBp, startBp + length, chr);
	}

	@Override
	public FileParser clone() {
		FileParser clone = new FastaParser(titleLength, rowLength, chr);
		
		clone.chunk = this.chunk;
			
		return clone;
	}
}