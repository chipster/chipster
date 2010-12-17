package fi.csc.microarray.gbrowser.index;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Getting genes information for indexing from a file
 * 
 */

public class GetGeneIndexData {

	private ChunkDataSource dataSource;

	public GetGeneIndexData(ChunkDataSource dataSource) {
		this.dataSource = dataSource;
	}

	public List<RegionContent> read() {

		try {
			byte[] fileChunk = dataSource.readAll();
			List<ColumnType> columns = Arrays.asList(new ColumnType[] { ColumnType.CHROMOSOME, ColumnType.BP_START, ColumnType.BP_END,
					ColumnType.DESCRIPTION });

			return dataSource.getFileParser().getAll(new Chunk(new String(fileChunk)), columns);

		} catch (FileNotFoundException e) {
			Session.getSession().getApplication().reportException(e);
			return null;
		} catch (IOException e) {
			Session.getSession().getApplication().reportException(e);
			return null;
		}
	}
}
