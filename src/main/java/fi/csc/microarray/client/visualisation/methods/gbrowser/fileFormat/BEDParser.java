package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.util.IOUtils;

/**
 * <p>Parser for BED file format.</p>
 * 
 * <p>BED files describe data lines that are displayed in an annotation track.
 * Information is stored using tab-separated values.<br/><br/>
 * 
 * Example:<br/>
 * <pre>
 * track name=pairedReads description="Clone Paired Reads" useScore=1
 * chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
 * chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
 * </pre></p>
 * 
 * @see http://genome.ucsc.edu/FAQ/FAQformat.html#format1
 * 
 * @author Petri Klemelä, Aleksi Kallio
 *
 */
public class BEDParser extends TsvParser {

	private static final String BED_HEADER_STRING = "track";

	public static List<ColumnDefinition> completeBedColumns = Arrays.asList(new ColumnDefinition[] { 
		new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING), 
		new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
		new ColumnDefinition(ColumnType.BP_END, Type.LONG),
		new ColumnDefinition(ColumnType.ID, Type.STRING),
		new ColumnDefinition(ColumnType.VALUE, Type.FLOAT),
		new ColumnDefinition(ColumnType.STRAND, Type.STRING),
		new ColumnDefinition(ColumnType.THICK_START, Type.STRING),
		new ColumnDefinition(ColumnType.THICK_END, Type.STRING),
		new ColumnDefinition(ColumnType.ITEM_RGB, Type.STRING),
		new ColumnDefinition(ColumnType.BLOCK_COUNT, Type.STRING),
		new ColumnDefinition(ColumnType.BLOCK_SIZES, Type.STRING),
		new ColumnDefinition(ColumnType.BLOCK_STARTS, Type.STRING),
	});

	public BEDParser() {
		super(new FileDefinition(completeBedColumns));
	}

	public BEDParser(FileDefinition fileDefinition) {
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
		
		if (col == ColumnType.BP_START || col == ColumnType.BP_END) {
			return (Long)obj + 1;
		}
		return obj;
	}
	
	public long getHeaderLength(String string) throws IOException {
		BufferedReader reader = new BufferedReader(new StringReader(string));
		return getHeaderLength(reader);
	}

	public long getHeaderLength(BufferedReader reader) throws IOException {
		String firstLine = reader.readLine();
		if (firstLine != null && firstLine.startsWith(BED_HEADER_STRING)) {
			return firstLine.length();
		} else {
			return 0;
		}
	}

	@Override
	public long getHeaderLength(File file) throws IOException {
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(new FileInputStream(file)));
			return getHeaderLength(in);
			
		} finally {
			IOUtils.closeIfPossible(in);
		}
	}
	
	@Override
	public BpCoordRegion getBpRegion(Chunk chunk) {
		BpCoordRegion reg = super.getBpRegion(chunk);
		return new BpCoordRegion(reg.start.bp + 1, reg.start.chr, reg.end.bp + 1, reg.end.chr);
	}
}