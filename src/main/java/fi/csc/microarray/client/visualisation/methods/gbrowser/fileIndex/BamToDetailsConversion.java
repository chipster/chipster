package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.BamDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.SingleThreadAreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;

/**
 * This conversion reads bam files with Picard and creates a RegionContent object for each read.
 * 
 * @author Aleksi Kallio, Petri Klemelä
 */
public class BamToDetailsConversion extends SingleThreadAreaRequestHandler {
	
	private static final int RESULT_CHUNK_SIZE = 100;

	private BamDataSource dataSource;

	private Region previousRequestedRegion;

	public BamToDetailsConversion(BamDataSource file, final GBrowser browser) {
	    
		super(null, null);
		
		this.dataSource = file;
	}
		
	@Override
	public void clean() {
		
		SamBamUtils.closeIfPossible(dataSource.getReader());
	}


	@Override
	protected void processAreaRequest(AreaRequest request) {
						
		super.processAreaRequest(request);
		
		if (request.getStatus().poison) {
			return;
		}
		
		if (request.getRequestedContents().contains(ColumnType.COVERAGE)) {
			
			request.getRequestedContents().remove(ColumnType.COVERAGE);
			
			request.getRequestedContents().addAll(Arrays.asList(new ColumnType[] {
			ColumnType.ID, 
			ColumnType.SEQUENCE,
			ColumnType.STRAND,
			ColumnType.QUALITY,
			ColumnType.CIGAR}));
			
			try {
				getDetails(request);
			} catch (IOException e) {
				e.printStackTrace();
			}		

		} else {
			
			try {
				getDetails(request);
			} catch (IOException e) {
				e.printStackTrace();
			}			
		}												
	}
		
	private void getDetails(AreaRequest request) throws IOException {
		
		fetchReads(request);
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
	public void fetchReads(AreaRequest request) {

		// Read the given region
		String chromosome = dataSource.getChromosomeNameUnnormaliser().unnormalise(request.start.chr);
		CloseableIterator<SAMRecord> iterator = dataSource.getReader().query(chromosome, request.start.bp.intValue(), request.end.bp.intValue(), false);

		// Produce results
		while (iterator.hasNext()) {

			List<RegionContent> responseList = new LinkedList<RegionContent>();

			// Split results into chunks
			for (int c = 0; c < RESULT_CHUNK_SIZE && iterator.hasNext(); c++) {
				SAMRecord record = iterator.next();

				// Region for this read
				Region recordRegion = new Region((long) record.getAlignmentStart(), (long) record.getAlignmentEnd(), request.start.chr);

				// Values for this read
				LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();

				RegionContent read = new RegionContent(recordRegion, values);

				if (request.getRequestedContents().contains(ColumnType.ID)) {
					values.put(ColumnType.ID, record.getReadName());
				}

				if (request.getRequestedContents().contains(ColumnType.STRAND)) {
					values.put(ColumnType.STRAND, record.getReadNegativeStrandFlag() ? Strand.REVERSE : Strand.FORWARD);
				}

				if (request.getRequestedContents().contains(ColumnType.QUALITY)) {
					values.put(ColumnType.QUALITY, record.getBaseQualityString());
				}

				if (request.getRequestedContents().contains(ColumnType.CIGAR)) {
					Cigar cigar = new Cigar(read, record.getCigar());
					values.put(ColumnType.CIGAR, cigar);
				}

				// TODO Deal with "=" and "N" in read string
				if (request.getRequestedContents().contains(ColumnType.SEQUENCE)) {
					String seq = record.getReadString();
					values.put(ColumnType.SEQUENCE, seq);
				}

				if (request.getRequestedContents().contains(ColumnType.MATE_POSITION)) {
					
					BpCoord mate = new BpCoord((Long)(long)record.getMateAlignmentStart(),
							new Chromosome(record.getMateReferenceName()));
					
					values.put(ColumnType.MATE_POSITION, mate);
				}
				
				/*
				 * NOTE! RegionContents created from the same read area has to be equal in methods equals, hash and compareTo. Primary types
				 * should be ok, but objects (including tables) has to be handled in those methods separately. Otherwise tracks keep adding
				 * the same reads to their read sets again and again.
				 */
				responseList.add(read);
			}

			// Send result			
			super.createAreaResult(new AreaResult(request.getStatus(), responseList));			
		}

		// We are done
		iterator.close();
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
