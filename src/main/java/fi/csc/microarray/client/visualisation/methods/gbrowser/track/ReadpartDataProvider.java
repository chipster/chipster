package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.util.Collection;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;

/**
 * Splices reads into spliced parts, returned as {@link ReadPart} objects.
 * 
 * @author Aleksi Kallio
 *
 */
public class ReadpartDataProvider implements AreaResultListener {

	private GBrowserView view;
	private AreaRequestHandler readData;
	private Collection<RegionContent> reads = new TreeSet<RegionContent>();
	private LinkedList<ReadPart> readParts = new LinkedList<ReadPart>(); 
	private LinkedList<ReadPart> readPartsF = new LinkedList<ReadPart>(); 
	private LinkedList<ReadPart> readPartsR = new LinkedList<ReadPart>();
	private boolean needsRefresh = false;

	public ReadpartDataProvider(GBrowserView view, AreaRequestHandler readData) {
		this.view = view;
		this.readData = readData;
		
		// start listening
		view.getQueueManager().addResultListener(readData, this);
	}

	@Override
	public void processAreaResult(AreaResult areaResult) {


		// Add this to queue of RegionContents to be processed
		synchronized (reads) {

			// Here identical region contents are removed (set semantics, no duplicates)
			// So it is essential that reads have their unique ID's.
			for (RegionContent read : areaResult.getContents()) {	  
				if (areaResult.getStatus().areaRequestHandler == readData && 
						read.values.containsKey(ColumnType.STRAND) &&
						read.values.containsKey(ColumnType.SEQUENCE) && 
						read.values.containsKey(ColumnType.CIGAR)) {
					
					this.reads.add(read);
					needsRefresh = true;
				}
			}
		}
		view.redraw();
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
			}
			throw new IllegalArgumentException("illegal strand: " + strand);
		}

	}

	private void refreshReadparts() {
		readParts.clear();
		readPartsF.clear();
		readPartsR.clear();
		Iterator<RegionContent> iter = reads.iterator();
		while (iter.hasNext()) {

			RegionContent read = iter.next();

			// Remove reads that are not in this view
			if (!view.requestIntersects(read.region)) {
				iter.remove();
				continue;
			}

			// Split read into continuous blocks (elements) by using the cigar
			List<ReadPart> visibleRegions = Cigar.splitElements(read);
			
			// Pool and sort read parts by strands
			for (ReadPart visibleRegion : visibleRegions) {
				// Skip read parts that are not in this view
				if (!visibleRegion.intersects(view.getBpRegion())) {
					continue;
				}
				
				readParts.add(visibleRegion); 
				
				if (read.values.get(ColumnType.STRAND) == Strand.FORWARD) {
					readPartsF.add(visibleRegion);
				} else if (read.values.get(ColumnType.STRAND) == Strand.REVERSE) {
					readPartsR.add(visibleRegion);
				}
			}
		}
		
		Collections.sort(readParts);
		Collections.sort(readPartsF);
		Collections.sort(readPartsR);
	}

}
