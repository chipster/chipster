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


public class VcfParser extends TsvParser {

	private static final String VCF_HEADER_START = "#";

	public static List<ColumnDefinition> completeVcfColumns = Arrays.asList(new ColumnDefinition[] { 
			new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING), 
			new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
	});

	private Long headerLength;

	public VcfParser() {
		super(new FileDefinition(completeVcfColumns));
	}

	public VcfParser(FileDefinition fileDefinition) {
		super(fileDefinition);
	}

	@Override
	public String getName() {
		return "Vcf";
	}

	@Override
	public long getDefaulChunkLength() {
		return 128;
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
		while (line != null && line.startsWith(VCF_HEADER_START)) {
			
			bytes += line.length() + 1; //plus one is the new line character
		} 

		return bytes;
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