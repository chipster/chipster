package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

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
	public RegionContent[] concise(String chunk) {

		long totalF = 0;
		long totalR = 0;
		
		long readLength = ((String)get(getFirstRow(chunk), ColumnType.SEQUENCE)).length();
						
		for (RegionContent rc : 
			getAll(chunk, Arrays.asList(new ColumnType[] { ColumnType.STRAND }))) {

			if ((Strand) rc.values.get(ColumnType.STRAND) == Strand.FORWARD) {
				totalF += readLength;

			} else {
				totalR += readLength;
			}
		}

		RegionContent[] result = new RegionContent[] {
				new RegionContent(getBpRegion(chunk), totalF / (float)getBpRegion(chunk).getLength()),
				new RegionContent(getBpRegion(chunk), totalR / (float)getBpRegion(chunk).getLength())
		};
		
		result[0].values.put(ColumnType.STRAND, Strand.FORWARD);
		result[1].values.put(ColumnType.STRAND, Strand.REVERSED);

		return result;
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