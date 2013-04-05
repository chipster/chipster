package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.CytobandDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRetrievalStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.QueueManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class ViewLimiter implements RegionListener {

	private QueueManager queueManager;
	private AreaRequestHandler cytobandDataSource;
	private BpCoord limit;
	
	private List<RegionListener> limitChangeListeners = new LinkedList<RegionListener>();

	/**
	 * @param queueManager
	 * @param cytobandRequestHandler
	 * @param view View to follow to notice chromosome changes
	 */
	public ViewLimiter(QueueManager queueManager, AreaRequestHandler cytobandRequestHandler, GBrowserView view) {
		this.queueManager = queueManager;
		this.cytobandDataSource = cytobandRequestHandler;

		queueManager.addResultListener(cytobandRequestHandler, new AreaResultListener() {

			@Override
			public void processAreaResult(AreaResult areaResult) {
				
				if (limit != null) {
					Long previousLimit = limit.bp;

					for (RegionContent regCont : areaResult.getContents()) {

						BpCoord value = regCont.region.end;

						if (value.chr.equals(limit.chr)) {
							if (value.bp > limit.bp) {
								limit.bp = value.bp;
							}
						} 
					}

					if (!previousLimit.equals(limit.bp)) {
						for (RegionListener listener : limitChangeListeners) {
							listener.regionChanged(new Region(0l, limit.bp, limit.chr));
						}
					}
				}
			}
		});

		view.addRegionListener(this);
	}


	public void regionChanged(Region bpRegion) {
		
		if (limit == null || !limit.chr.equals(bpRegion.start.chr)) {

			limit = new BpCoord(0l, bpRegion.start.chr);

			queueManager.addAreaRequest(
					cytobandDataSource, new AreaRequest(new Region(0l, Long.MAX_VALUE, bpRegion.start.chr), 
							new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {ColumnType.VALUE })), 
							new DataRetrievalStatus()), false);
		}
	}

	public BpCoord getLimit() {
		if (limit != null) {
			return limit.clone();
		} else {
			return null;
		}
	}
	
	public void addLimitChangeListener(RegionListener listener) {
		limitChangeListeners.add(listener);
	}
	
	public void removeLimitChangeListeners(RegionListener listener) {
		limitChangeListeners.remove(listener);
	}
}
