package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserSettings.CoverageType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Base;
import fi.csc.microarray.util.BamUtils;

/**
 * This conversion class uses Picard to read Bam files and calculates a coverage.
 * 
 * @author Petri Klemelä
 *
 */
public class BamToCoverageConversion extends DataThread {
	
	private BamDataSource dataSource;

	private CoverageType coverageType;
	
	public BamToCoverageConversion(BamDataSource file, CoverageType coverageType, final GBrowser browser) {
	    
		super(browser, file);
		
		this.dataSource = file;
		this.coverageType = coverageType;
	}
		
	@Override
	public void clean() {
		
		dataSource.close();
	}


	@Override
	protected void processDataRequest(DataRequest request) throws InterruptedException {							
		
		long step = 10_000;
		
		// Divide visible region into subregions and iterate over them
		for (long pos = request.start.bp; pos < request.end.bp; pos += step ) {

			BpCoord from = new BpCoord(pos, request.start.chr);
			BpCoord to = new BpCoord(Math.min(pos + step, request.end.bp), request.start.chr);

			calculateCoverage(request, from, to);									
		}
	}

	private void calculateCoverage(DataRequest request, BpCoord from, BpCoord to) throws InterruptedException {	
				
		//query data for full average bins, because merging them later would be difficult
		long start = CoverageTool.getBin(from.bp);		
		long end = CoverageTool.getBin(to.bp) + CoverageTool.BIN_SIZE - 1;

		//if end coordinate is smaller than 1 Picard returns a full chromosome and we'll run out of memory
		if (end < 1) {
			end = 1;
		}
				
		CloseableIterator<SAMRecord> iterator = dataSource.query(from.chr, (int)start, (int)end);		
		
		BaseStorage forwardBaseStorage = new BaseStorage();
		BaseStorage reverseBaseStorage = new BaseStorage();
		
		while (iterator.hasNext()) {
			
			SAMRecord record = iterator.next();
			
			LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();

			Region recordRegion = new Region((long) record.getAlignmentStart(), (long) record.getAlignmentEnd(), request.start.chr);
			
			Feature read = new Feature(recordRegion, values);

			values.put(DataType.ID, record.getReadName());
			
			values.put(DataType.STRAND, BamUtils.getStrand(record, coverageType));
			
			Cigar cigar = new Cigar(read, record.getCigar());
			values.put(DataType.CIGAR, cigar);
			
			String seq = record.getReadString();
			values.put(DataType.SEQUENCE, seq);
			
			
			// Split read into continuous blocks (elements) by using the cigar
			List<ReadPart> parts = Cigar.splitElements(read);
			
			for (ReadPart part : parts) {				 
				
				if (read.values.get(DataType.STRAND) == Strand.FORWARD) {
					forwardBaseStorage.addNucleotideCounts(part);
				} else if (read.values.get(DataType.STRAND) == Strand.REVERSE) {
					reverseBaseStorage.addNucleotideCounts(part);			
				}
			}						
		}
				
		// We are done
		iterator.close();
				
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
