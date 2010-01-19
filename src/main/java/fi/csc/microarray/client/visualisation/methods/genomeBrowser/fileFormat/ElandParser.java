package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public class ElandParser extends ConstantRowLengthParser{

	public ElandParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {

								new ColumnDefinition(ColumnType.ID, Type.STRING, 32),
								new ColumnDefinition(ColumnType.SEQUENCE, Type.STRING, 64),
								new ColumnDefinition(ColumnType.QUALITY, Type.STRING, 8),
								new ColumnDefinition(ColumnType.SKIP, Type.STRING, 8),
								new ColumnDefinition(ColumnType.SKIP, Type.STRING, 8),
								new ColumnDefinition(ColumnType.SKIP, Type.STRING, 8),
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING, 16),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG, 16),							
								new ColumnDefinition(ColumnType.STRAND, Type.STRING, 2),
								new ColumnDefinition(ColumnType.SKIP, Type.NEWLINE, 1),

						})));
	}

	@Override
	public int getChunkMaxByteLength() {
		return (int)getRowByteLength() * 32;
	}

	@Override
	public long getFilePosition(long readIndex) {
		return readIndex * getRowByteLength();
	}

	@Override
	public long getRowIndex(long filePosition) {
		return filePosition / getRowByteLength();
	}

	@Override
	public RegionContent[] concise(BpCoordRegion bpRegion) {

		ColumnType start = ColumnType.BP_START;
		ColumnType seq = ColumnType.SEQUENCE;

		long totalF = 0;
		long totalR = 0;
		
		int i;
		long minBp = -1;
		long maxBp = -1;

		long rowCount = getChunkRowCount();

		long length = 0;
		
		for (i = 0; i < rowCount; i++){
			
			if(i == 0){				
				length = ((String)get(i + chunk.rowIndex, seq)).length();
			}

			long startBp = (Long)get(i + chunk.rowIndex, start);

			if((Strand)get(i + chunk.rowIndex, ColumnType.STRAND) == Strand.FORWARD) {
				
				totalF += length;
				
			} else {
				totalR += length;
			}
			
			if(i == 0){
				minBp = startBp;
			}
			if(i == rowCount - 1){
				maxBp = startBp;			
			}
		}
		
		RegionContent[] result = new RegionContent[] {
			new RegionContent(bpRegion, totalF / (float)(maxBp - minBp)),
			new RegionContent(bpRegion, totalR / (float)(maxBp - minBp))
		};		
		
		result[0].values.put(ColumnType.STRAND, Strand.FORWARD);
		result[1].values.put(ColumnType.STRAND, Strand.REVERSED);
		
		return result;
	}

	@Override
	public BpCoordRegion getBpRegion(long rowIndex) {

		long startBp = (Long)get(rowIndex, ColumnType.BP_START);
		long length = ((String)get(rowIndex, ColumnType.SEQUENCE)).trim().length();
		Chromosome chr = (Chromosome)get(rowIndex, ColumnType.CHROMOSOME);

		return new BpCoordRegion(startBp, startBp + length, chr);
	}

	@Override
	public FileParser clone() {
		FileParser clone = new ElandParser();
		
		clone.chunk = this.chunk;
		
		return clone;
	}
	
	@Override
	public Object get(long rowIndex, ColumnType col) {
		
		Object obj = super.get(rowIndex, col);
		
		if(col == ColumnType.CHROMOSOME) {			
			return new Chromosome(((Chromosome)obj).toString().replace(".fa", ""));
		}
		return obj;
	}
}