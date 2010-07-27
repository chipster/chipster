package fi.csc.microarray.gbrowser.index;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.GeneParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Getting genes information for indexing from a file
 */

public class GetGeneIndexData {
	
	private static final File URL_ROOT;

	static {
		URL_ROOT = new File("/home/zukauska/chipster-share/ngs/annotations");
	}
	
	public GetGeneIndexData (){
		
	}
	
	public List<RegionContent> read(){
        
        try {
        	File geneFile = new File(URL_ROOT+"/Homo_sapiens.NCBI36.54_genes.tsv");
			ChunkDataSource data = new ChunkDataSource(
					geneFile,
					new GeneParser());
			byte[] fileChunk = new byte[(int)geneFile.length()];
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
