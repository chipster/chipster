package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.SortedMap;
import java.util.TreeMap;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class ReadpartDataProvider implements AreaResultListener {

	private View view;
	private DataSource readData;
	private Collection<RegionContent> reads = new TreeSet<RegionContent>();
	private SortedMap<Long, ReadPart> readParts = new TreeMap<Long, ReadPart>(); 
	private SortedMap<Long, ReadPart> readPartsF = new TreeMap<Long, ReadPart>(); 
	private SortedMap<Long, ReadPart> readPartsR = new TreeMap<Long, ReadPart>(); 

	public ReadpartDataProvider(View view, DataSource readData, Class<? extends AreaRequestHandler> readDataHandler) {
		this.view = view;
		this.readData = readData;
		
		// start listening
		view.getQueueManager().createQueue(readData, readDataHandler);
		view.getQueueManager().addResultListener(readData, this);
	}

	@Override
	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		// Check that areaResult has false concised status and correct strand
		if (areaResult.status.file == readData && areaResult.status.concise == false) {
			// Add this to queue of RegionContents to be processed
			this.reads.add(areaResult.content);
			refreshReadparts();
			view.redraw();
		}
	}

	private void refreshReadparts() {
		// TODO Auto-generated method stub
		 //&& areaResult.content.values.get(ColumnType.STRAND) == getStrand()
		
		readParts.clear();
		readPartsF.clear();
		readPartsR.clear();
		Iterator<RegionContent> iter = reads.iterator();
		while (iter.hasNext()) {

			RegionContent read = iter.next();

			// Remove reads that are not in this view
			if (!read.region.intersects(view.getBpRegion())) {
				iter.remove();
				continue;
			}

			// Split read into continuous blocks (elements) by using the cigar
			List<ReadPart> visibleRegions = Cigar.splitVisibleElements(read);
			
			// Pool and sort read parts by strands
			for (ReadPart visibleRegion : visibleRegions) {
				readParts.put(visibleRegion.start.bp, visibleRegion); 
				
				if (read.values.get(ColumnType.STRAND) == Strand.FORWARD) {
					readPartsF.put(visibleRegion.start.bp, visibleRegion);
				} else if (read.values.get(ColumnType.STRAND) == Strand.REVERSED) {
					readPartsR.put(visibleRegion.start.bp, visibleRegion);
				}

			}
		}
	}

	public SortedMap<Long, ReadPart> getReadparts(Strand strand) {
		switch (strand) {
			case BOTH:
				return readParts;
			case FORWARD:
				return readPartsF;
			case REVERSED:
				return readPartsR;
		}
		throw new IllegalArgumentException("illegal strand: " + strand);
	}

}
