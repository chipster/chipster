package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;

/**
 * Splices reads into spliced parts, returned as {@link ReadPart} objects.
 * 
 * @author Aleksi Kallio
 *
 */
public class ReadpartDataProvider {

	private Collection<Feature> reads = new TreeSet<Feature>();
	private LinkedList<ReadPart> readParts = new LinkedList<ReadPart>(); 
	private LinkedList<ReadPart> readPartsF = new LinkedList<ReadPart>(); 
	private LinkedList<ReadPart> readPartsR = new LinkedList<ReadPart>();
	private boolean needsRefresh = false;
	
	public void addReads(LinkedList<Feature> reads) {


		// Add this to queue of RegionContents to be processed
		synchronized (reads) {

			// Here identical region contents are removed (set semantics, no duplicates)
			// So it is essential that reads have their unique ID's.
			for (Feature read : reads) {	  
				if (
						read.values.containsKey(DataType.STRAND) &&
						read.values.containsKey(DataType.SEQUENCE) && 
						read.values.containsKey(DataType.CIGAR)) {
					
					this.reads.add(read);
					needsRefresh = true;
				}
			}
		}
	}

	public Iterable<ReadPart> getReadparts(Strand strand) {

		synchronized (reads) {

			if (needsRefresh) {
				refreshReadparts();
				needsRefresh = false;
			}

			switch (strand) {
			case BOTH:
				return readParts;
			case FORWARD:
				return readPartsF;
			case REVERSE:
				return readPartsR;
			default:
				throw new IllegalArgumentException("illegal strand: " + strand);
			}
		}
	}

	private void refreshReadparts() {
		readParts.clear();
		readPartsF.clear();
		readPartsR.clear();
		Iterator<Feature> iter = reads.iterator();
		while (iter.hasNext()) {

			Feature read = iter.next();

			// Split read into continuous blocks (elements) by using the cigar
			List<ReadPart> visibleRegions = Cigar.splitElements(read);
			
			// Pool and sort read parts by strands
			for (ReadPart visibleRegion : visibleRegions) {
				// Skip read parts that are not in this view				
				readParts.add(visibleRegion); 
				
				if (read.values.get(DataType.STRAND) == Strand.FORWARD) {
					readPartsF.add(visibleRegion);
				} else if (read.values.get(DataType.STRAND) == Strand.REVERSE) {
					readPartsR.add(visibleRegion);
				}
			}
		}
		
		Collections.sort(readParts);
		Collections.sort(readPartsF);
		Collections.sort(readPartsR);
	}
}
