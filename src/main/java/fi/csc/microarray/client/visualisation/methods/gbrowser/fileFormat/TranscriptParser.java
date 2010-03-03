package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class TranscriptParser extends ConstantRowLengthParser{

	public TranscriptParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {

								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG, 2),
								new ColumnDefinition(ColumnType.PARENT_BP_START, Type.LONG, 16),
								new ColumnDefinition(ColumnType.PARENT_BP_END, Type.LONG, 16),
								new ColumnDefinition(ColumnType.STRAND, Type.STRING, 2),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG, 16),
								new ColumnDefinition(ColumnType.BP_END, Type.LONG, 16),
								new ColumnDefinition(ColumnType.DESCRIPTION, Type.STRING, 32),
								new ColumnDefinition(ColumnType.VALUE, Type.STRING, 16),							
								new ColumnDefinition(ColumnType.PARENT_ID, Type.STRING, 32),
								
								//Should be Type.LONG, but there is some quotes in current annotation file
								new ColumnDefinition(ColumnType.PARENT_PART, Type.STRING, 4),
								new ColumnDefinition(ColumnType.SKIP, Type.NEWLINE, 1)

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

		long totalF = 0;
		long totalR = 0;
		
		int i;
		long minBp = -1;
		long maxBp = -1;

		long rowCount = getChunkRowCount();
		
		for (i = 0; i < rowCount; i++){			

			long startBp = (Long)get(i + chunk.rowIndex, ColumnType.BP_START);
			long endBp = (Long)get(i + chunk.rowIndex, ColumnType.BP_END);
			long length = endBp - startBp + 1;

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
		long endBp = (Long)get(rowIndex, ColumnType.BP_END);

		Chromosome chr = (Chromosome)get(rowIndex, ColumnType.CHROMOSOME);

		return new BpCoordRegion(startBp, endBp, chr);
	}

	@Override
	public FileParser clone() {
		FileParser clone = new TranscriptParser();
		
		clone.chunk = this.chunk;
		
		return clone;
	}

	@Override
	public String getName() {
		return "Chipster transcript file";
	}
}