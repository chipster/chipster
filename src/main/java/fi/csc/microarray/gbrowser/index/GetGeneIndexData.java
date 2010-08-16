package fi.csc.microarray.gbrowser.index;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.GeneParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Getting genes information for indexing from a file
 */

public class GetGeneIndexData {
	
	private GenomeBrowser genomeBrowser;
	
	public GetGeneIndexData(GenomeBrowser genomeBrowser) {
		this.genomeBrowser = genomeBrowser;
	}
	
	public List<RegionContent> read(){
        
        try {
			ChunkDataSource data = genomeBrowser.createAnnotationDataSource(
			        "Homo_sapiens.NCBI36.54_genes.tsv", new GeneParser());
			byte[] fileChunk = new byte[(int)data.getSize()];
			data.read(0, fileChunk);
			List<ColumnType> columns = Arrays.asList(new ColumnType[] {
					ColumnType.CHROMOSOME,
					ColumnType.BP_START,
					ColumnType.BP_END,
					ColumnType.DESCRIPTION});
			
			return  data.getFileParser().getAll(new Chunk(new String(fileChunk)), columns);
			
		} catch (FileNotFoundException e) {
			return null;
		} catch (IOException e) {
			return null;
		}
	}
}
