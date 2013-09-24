package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadpartDataProvider;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Base;

/**
 * This conversion class uses Picard to read Bam files and calculates a coverage.
 * 
 * @author Petri Klemelä
 *
 */
public class BamToCoverageConversion extends DataThread {
	
	private BamDataSource dataSource;
	
	public BamToCoverageConversion(BamDataSource file, final GBrowser browser) {
	    
		super(browser, file);
		
		this.dataSource = file;
	}
		
	@Override
	public void clean() {
		
		dataSource.close();
	}


	@Override
	protected void processDataRequest(DataRequest request) {							
		
		long step = 1*10000;
		
		// Divide visible region into subregions and iterate over them
		for (long pos = request.start.bp; pos < request.end.bp; pos += step ) {

			BpCoord from = new BpCoord(pos, request.start.chr);
			BpCoord to = new BpCoord(Math.min(pos + step, request.end.bp), request.start.chr);

			calculateCoverage(request, from, to);									
		}
	}

	private void calculateCoverage(DataRequest request, BpCoord from, BpCoord to) {	
		
		//query data for full average bins, because merging them later would be difficult
		long start = CoverageTool.getBin(from.bp);
		long end = CoverageTool.getBin(to.bp) + CoverageTool.BIN_SIZE - 1;
				
		CloseableIterator<SAMRecord> iterator = dataSource.query(from.chr, (int)start, (int)end);
		
		LinkedList<Feature> reads = new LinkedList<Feature>();
		
		while (iterator.hasNext()) {
			
			SAMRecord record = iterator.next();
			
			LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();

			Region recordRegion = new Region((long) record.getAlignmentStart(), (long) record.getAlignmentEnd(), request.start.chr);
			Feature read = new Feature(recordRegion, values);

			values.put(DataType.ID, record.getReadName());
			
			values.put(DataType.STRAND, record.getReadNegativeStrandFlag() ? Strand.REVERSE : Strand.FORWARD);
			
			Cigar cigar = new Cigar(read, record.getCigar());
			values.put(DataType.CIGAR, cigar);
			
			String seq = record.getReadString();
			values.put(DataType.SEQUENCE, seq);
			
			reads.add(read);
		}				
				
		// We are done
		iterator.close();
		
		ReadpartDataProvider readpartProvider = new ReadpartDataProvider();
		readpartProvider.addReads(reads);		
		BaseStorage forwardBaseStorage = new BaseStorage();
		BaseStorage reverseBaseStorage = new BaseStorage();
		
		forwardBaseStorage.getNucleotideCounts(readpartProvider.getReadparts(Strand.FORWARD), null, null);
		reverseBaseStorage.getNucleotideCounts(readpartProvider.getReadparts(Strand.REVERSE), null, null);
		
		/* Reads that overlap query regions create nucleotide counts outside the query region.
		 * Remove those extra nucleotide counts, because they don't contain all reads of those regions and would show
		 * wrong information. 
		 */
		Region filterRegion = new Region(start, end, from.chr);
		forwardBaseStorage.filter(filterRegion);
		reverseBaseStorage.filter(filterRegion);

		// Send result		
		LinkedList<Feature> resultList = new LinkedList<Feature>();					
		
		createResultList(from, forwardBaseStorage, resultList, Strand.FORWARD);
		createResultList(from, reverseBaseStorage, resultList, Strand.REVERSE);		
		
		LinkedList<Feature> averageCoverage = CoverageTool.average(resultList, from.chr);
		
		if (request.getRequestedContents().contains(DataType.COVERAGE)) {		
			super.createDataResult(new DataResult(request, resultList));
		}
		
		if (request.getRequestedContents().contains(DataType.COVERAGE_AVERAGE)) {
			super.createDataResult(new DataResult(request, averageCoverage));
		}
	}

	private void createResultList(BpCoord from, BaseStorage baseStorage,
			LinkedList<Feature> resultList, Strand strand) {
		
		Iterator<Base> baseIter = baseStorage.iterator();
		
		while(baseIter.hasNext()) {
			Base base = baseIter.next();
			Region region = new Region(base.getBpLocation(), base.getBpLocation(), from.chr);
			
			LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();
			values.put(DataType.VALUE, base);
			values.put(DataType.STRAND, strand);
			
			resultList.add(new Feature(region, values));
		}
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
