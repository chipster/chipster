package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.VcfParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Very limited tool for handling vcf files similarly to how RegionOperations handles bed files.
 * 
 * @author Petri Klemelä
 *
 */
public class PositionOperations {

	/**
	 * Parses regions from a BED text formatted input file.
	 * 
	 * @param input BED file
	 * @return regions and their extra data
	 * @throws URISyntaxException 
	 */
	public List<RegionContent> loadFile(File input) throws FileNotFoundException, IOException, URISyntaxException {
		ChunkDataSource dataSource = new ChunkDataSource(input.toURI().toURL(), new VcfParser(), ChunkTreeHandlerThread.class);
		byte[] fileChunk = dataSource.readAll();
		return parseString(new String(fileChunk));
	}

	/**
	 * Parses regions from a BED text formatted String.
	 * 
	 * @param string BED string
	 * @return  regions and their extra data
	 */
	public List<RegionContent> parseString(String string) throws FileNotFoundException, IOException {
		
		// Process track name, if exists
		VcfParser parser = new VcfParser();
		int headerLength = (int)parser.getHeaderLength(string);
		string = string.substring(headerLength > 0 ? (headerLength + 1) : 0);
		
		// Count fields and create list of what extra types we need
		int fieldCount = string.split("\n")[0].split("\t").length;
		if (fieldCount < 2) {
			throw new IllegalArgumentException("VCF must have at least chromosome and position fields");
		}
		
		LinkedList<ColumnType> extraTypes = new LinkedList<ColumnType>();
		
// There aren't other fields defined in VcfParser
//		for (int i = 2; i < fieldCount; i++) {
//			extraTypes.add(VcfParser.completeVcfColumns.get(i).content);
//		}
		
		// Parse it
		return parser.getAll(new Chunk(string), extraTypes);
	}
}
