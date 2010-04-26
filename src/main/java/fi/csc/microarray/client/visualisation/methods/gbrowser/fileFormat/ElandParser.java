package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class ElandParser extends TsvParser {

	public ElandParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {
								new ColumnDefinition(ColumnType.ID, Type.STRING),
								new ColumnDefinition(ColumnType.SEQUENCE, Type.STRING),
								new ColumnDefinition(ColumnType.QUALITY, Type.STRING),
								new ColumnDefinition(ColumnType.SKIP, Type.STRING),
								new ColumnDefinition(ColumnType.SKIP, Type.STRING),
								new ColumnDefinition(ColumnType.SKIP, Type.STRING),
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG),							
								new ColumnDefinition(ColumnType.STRAND, Type.STRING)
						})));
	}
		
	@Override
	public RegionContent[] concise(BpCoordRegion nodeBpRegion) {

		long totalF = 0;
		long totalR = 0;
		
		long length = ((String)get(getFirstRow(), ColumnType.SEQUENCE)).length();
						
		for (RegionContent rc : 
			getAll(Arrays.asList(new ColumnType[] { ColumnType.STRAND }))) {

			if ((Strand) rc.values.get(ColumnType.STRAND) == Strand.FORWARD) {
				totalF += length;

			} else {
				totalR += length;
			}
		}

		RegionContent[] result = new RegionContent[] { 
				new RegionContent(nodeBpRegion, totalF / (float) nodeBpRegion.getLength()), 
				new RegionContent(nodeBpRegion, totalR / (float) nodeBpRegion.getLength()) };
		
		result[0].values.put(ColumnType.STRAND, Strand.FORWARD);
		result[1].values.put(ColumnType.STRAND, Strand.REVERSED);

		return result;
	}

	@Override
	public FileParser clone() {
		FileParser clone = new ElandParser();

		clone.chunk = this.chunk;

		return clone;
	}

	@Override
	public Object get(String[] cols, ColumnType col) {

		Object obj = super.get(cols, col);

		if (col == ColumnType.CHROMOSOME) {
			return new Chromosome(((Chromosome) obj).toString().replace(".fa", ""));
		}
		
		return obj;
	}

	@Override
	public String getName() {
		return "eland";
	}
}