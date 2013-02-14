package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.Arrays;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.util.IOUtils;

public class GtfParser extends TsvParser {

	private static final String GTF_HEADER_START = "#";

	public static List<ColumnDefinition> columns = Arrays.asList(new ColumnDefinition[] {
			new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING), //seqname
			new ColumnDefinition(ColumnType.SKIP, Type.STRING), //source
			new ColumnDefinition(ColumnType.VALUE, Type.STRING), //feature
			new ColumnDefinition(ColumnType.BP_START, Type.LONG), //start
			new ColumnDefinition(ColumnType.BP_END, Type.LONG), //end
			new ColumnDefinition(ColumnType.QUALITY, Type.STRING), //score
			new ColumnDefinition(ColumnType.STRAND, Type.STRING), //
			new ColumnDefinition(ColumnType.SKIP, Type.STRING), //frame
			new ColumnDefinition(ColumnType.METADATA, Type.STRING), //attributes
	});

	private Long headerLength;

	public GtfParser() {
		super(new FileDefinition(columns));
	}

	public GtfParser(FileDefinition fileDefinition) {
		super(fileDefinition);
	}

	@Override
	public String getName() {
		return "Gtf";
	}

	@Override
	public long getDefaulChunkLength() {
		return 256;
	}

	@Override
	public RegionContent[] concise(Chunk chunk) {
		return new RegionContent[] {};
	}

	public long getHeaderLength(String string) throws IOException {
		BufferedReader reader = new BufferedReader(new StringReader(string));
		return getHeaderLength(reader);
	}

	public long getHeaderLength(BufferedReader reader) throws IOException {
		String line = reader.readLine();
		long bytes = 0;
		while (line != null && line.startsWith(GTF_HEADER_START)) {
			
			bytes += line.length() + 1; //plus one for the new line character
			
			line = reader.readLine();
		} 

		return bytes - 1; //TODO find out why we have to subtract one to avoid losing first line of the real content
	}

	@Override
	public long getHeaderLength(File file) throws IOException {

		if (this.headerLength == null) {
			BufferedReader in = null;
			try {
				in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
				this.headerLength = getHeaderLength(in);

			} finally {
				IOUtils.closeIfPossible(in);
			}
		}
		return headerLength;
	}
}