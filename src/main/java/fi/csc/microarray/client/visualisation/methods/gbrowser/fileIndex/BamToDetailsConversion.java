package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

/**
 * This conversion reads bam files with Picard and creates a RegionContent object for each read.
 * 
 * @author Aleksi Kallio, Petri Klemelä
 */
public class BamToDetailsConversion extends DataThread {
	
	private static final int RESULT_CHUNK_SIZE = 100;

	private BamDataSource dataSource;

	public BamToDetailsConversion(BamDataSource file, final GBrowser browser) {
	    
		super(browser, file);
		
		this.dataSource = file;
	}
		
	@Override
	public void clean() {
		
		dataSource.close();
	}
		
	/**
	 * Find reads in a given range.
	 * 
	 * <p>
	 * TODO add cigar to the list of returned values
	 * <p>
	 * TODO add pair information to the list of returned values
	 * 
	 * @param request
	 * @return
	 */
	@Override
	protected void processDataRequest(DataRequest request) {

		// Read the given region
		CloseableIterator<SAMRecord> iterator = dataSource.query(request.start.chr, request.start.bp.intValue(), request.end.bp.intValue());		

		// Produce results
		while (iterator.hasNext()) {

			List<Feature> responseList = new LinkedList<Feature>();

			// Split results into chunks
			for (int c = 0; c < RESULT_CHUNK_SIZE && iterator.hasNext(); c++) {
				SAMRecord record = iterator.next();

				// Region for this read
				Region recordRegion = new Region((long) record.getAlignmentStart(), (long) record.getAlignmentEnd(), request.start.chr);

				// Values for this read
				LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();

				Feature read = new Feature(recordRegion, values);

				if (request.getRequestedContents().contains(DataType.ID)) {
					values.put(DataType.ID, record.getReadName());
				}

				if (request.getRequestedContents().contains(DataType.STRAND)) {
					values.put(DataType.STRAND, record.getReadNegativeStrandFlag() ? Strand.REVERSE : Strand.FORWARD);
				}

				if (request.getRequestedContents().contains(DataType.QUALITY)) {
					values.put(DataType.QUALITY, record.getBaseQualityString());
				}

				if (request.getRequestedContents().contains(DataType.CIGAR)) {
					Cigar cigar = new Cigar(read, record.getCigar());
					values.put(DataType.CIGAR, cigar);
				}

				// TODO Deal with "=" and "N" in read string
				if (request.getRequestedContents().contains(DataType.SEQUENCE)) {
					String seq = record.getReadString();
					values.put(DataType.SEQUENCE, seq);
				}

				if (request.getRequestedContents().contains(DataType.MATE_POSITION)) {
					
					BpCoord mate = new BpCoord((Long)(long)record.getMateAlignmentStart(),
							new Chromosome(record.getMateReferenceName()));
					
					values.put(DataType.MATE_POSITION, mate);
				}
				
				/*
				 * NOTE! RegionContents created from the same read area has to be equal in methods equals, hash and compareTo. Primary types
				 * should be ok, but objects (including tables) has to be handled in those methods separately. Otherwise tracks keep adding
				 * the same reads to their read sets again and again.
				 */
				responseList.add(read);
			}

			// Send result			
			super.createDataResult(new DataResult(request.getStatus(), responseList));			
		}

		// We are done
		iterator.close();
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
