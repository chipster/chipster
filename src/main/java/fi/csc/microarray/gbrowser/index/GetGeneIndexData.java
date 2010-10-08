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
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Getting genes information for indexing from a file
 * 
 * @author zukauska
 */

public class GetGeneIndexData {
	
	private GenomeBrowser genomeBrowser;
	private Genome selectedGenome;
	
	public GetGeneIndexData(GenomeBrowser genomeBrowser, Genome selectedGenome) {
		this.genomeBrowser = genomeBrowser;
		this.selectedGenome = selectedGenome;
	}
	
	public List<RegionContent> read(){
        
        try {
			ChunkDataSource data = genomeBrowser.createAnnotationDataSource(
					selectedGenome.species.replace(" ", "_") + "." + 
					selectedGenome.version + "_genes.tsv", new GeneParser());
			byte[] fileChunk = data.readAll();
			List<ColumnType> columns = Arrays.asList(new ColumnType[] {
					ColumnType.CHROMOSOME,
					ColumnType.BP_START,
					ColumnType.BP_END,
					ColumnType.DESCRIPTION});
			
			return  data.getFileParser().getAll(new Chunk(new String(fileChunk)), columns);
			
		} catch (FileNotFoundException e) {
			genomeBrowser.getClientApplication().reportException(e);
			return null;
		} catch (IOException e) {
			genomeBrowser.getClientApplication().reportException(e);
			return null;
		}
	}
}
