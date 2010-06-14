package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeqBlockTrack;

/**
 * Parser for a BED-like format that describes annotations. Values
 * are tab-separated.
 * 
 * @see BEDParser
 *
 */
public class BEDReadParser extends TsvParser {

	public BEDReadParser() {
		super(new FileDefinition(Arrays.asList(
				new ColumnDefinition[] { 
						new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING), 
						new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
						new ColumnDefinition(ColumnType.BP_END, Type.LONG), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING),
						new ColumnDefinition(ColumnType.STRAND, Type.STRING),
						new ColumnDefinition(ColumnType.QUALITY, Type.STRING),
						}
				)));
	}

	public BEDReadParser(FileDefinition fileDefinition) {
		super(fileDefinition);
	}

	@Override
	public String getName() {
		return "Chipster peaks";
	}

	@Override
	public long getDefaulChunkLength() {
		return 128;
	}

	@Override
	public RegionContent[] concise(String chunk) {
		return new RegionContent[] {};
	}

	@Override
	public Object get(String[] cols, ColumnType col) {

		if (col == ColumnType.SEQUENCE) {
			Long start = (Long)super.get(cols, ColumnType.BP_START);
			Long end = (Long)super.get(cols, ColumnType.BP_END);
			
			return SeqBlockTrack.DUMMY_SEQUENCE.substring(0, (int)(end-start));
			
		} else {

			Object obj = super.get(cols, col);

			if (col == ColumnType.CHROMOSOME) {
				return new Chromosome(((Chromosome) obj).toString().replace(".fa", ""));
			}

			return obj;
		}
	}
}