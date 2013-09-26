package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.IOException;
import java.net.URISyntaxException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Interpretation.TrackType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * @author klemela
 */
public class TsvLineParser extends AbstractTsvLineParser {		 

	private String[] header;
	private int chrColumn;
	private int startColumn;
	private int endColumn;
	private String headerStart;
	

	public TsvLineParser(DataUrl data, TrackType trackType) throws IOException, URISyntaxException {
		this(data);
		
		if (TrackType.TSV == trackType) {
			setChrColumn(0);			
		} else if (TrackType.TSV_WITH_ROW_ID == trackType) {
			setChrColumn(1);
		} else {
			throw new IllegalArgumentException("Unsupported trackType");
		}
	}
	
	public TsvLineParser(DataUrl data, int chrColumn) throws IOException, URISyntaxException {
		
		this(data);	
		setChrColumn(chrColumn);				
	}
	
	private TsvLineParser(DataUrl data) throws IOException, URISyntaxException {
		
		LineDataSource dataSource = new LineDataSource(data);
		headerStart = dataSource.readLine();
		String contentRow = dataSource.readLine();
		
		String[] splittedHeader = headerStart.split("\t");
		String[] splittedContent = contentRow.split("\t");
		
		if (splittedHeader.length == splittedContent.length - 1) {
			//Some R results don't have column title for the first column. (similar to generic Chipster implementation in TableColumnProvider)
			String[] headerWithId = new String[contentRow.length()];
			headerWithId[0] = "";
			System.arraycopy(splittedHeader, 0, headerWithId, 1, splittedHeader.length);
			
			this.header = headerWithId;  
		} else {
			this.header = splittedHeader;
		}
		
	}
	
	private void setChrColumn(int chrColumn) {
		this.chrColumn = chrColumn;
		this.startColumn = chrColumn + 1;
		this.endColumn = chrColumn + 2;
	}

	@Override
	public Region getRegion() {
		
		if (isContentLine()) {
			
			long start = getLong(startColumn);
			long end = getLong(endColumn);
			
			Chromosome chr = new Chromosome(getString(chrColumn));
			return new Region(start, end, chr);
			
		} else {
			//This is header line
			return null;
		}
	}		
	
	@Override
	public String getHeaderStart() {
		return headerStart;
	}	

	@Override
	public FileLine getFileLine() {
		TsvLine line = new TsvLine();
		
		line.setRegion(getRegion());
		line.setHeaders(header);
		line.setValues(values);
		
		return line;
	}
}
