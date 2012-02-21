package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Extends {@link TsvParser} to support files that have header line as their first line.
 *  
 * @author Aleksi Kallio
 *
 */
public class HeaderTsvParser extends TsvParser {

	public HeaderTsvParser() {
		super(new FileDefinition(Arrays.asList(
				new ColumnDefinition[] { 
						new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING), 
						new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
						new ColumnDefinition(ColumnType.BP_END, Type.LONG), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING), 
				})));
	}

	public HeaderTsvParser(FileDefinition fileDefinition) {
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
	public RegionContent[] concise(Chunk chunk) {
		return new RegionContent[] {};
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
	public long getHeaderLength(File file) throws IOException {
		BufferedReader in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
		String firstLine = in.readLine();
		return (firstLine + "\n").length();
	}
}