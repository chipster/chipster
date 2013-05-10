package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

public class ViewLimiter implements RegionListener {

	private QueueManager queueManager;
	private DataThread cytobandDataSource;
	private BpCoord limit;
	
	private List<RegionListener> limitChangeListeners = new LinkedList<RegionListener>();

	/**
	 * @param queueManager
	 * @param cytobandRequestHandler
	 * @param view View to follow to notice chromosome changes
	 */
	public ViewLimiter(QueueManager queueManager, DataThread cytobandRequestHandler, GBrowserView view) {
		this.queueManager = queueManager;
		this.cytobandDataSource = cytobandRequestHandler;

		queueManager.addDataResultListener(cytobandRequestHandler, new DataResultListener() {

			@Override
			public void processDataResult(DataResult dataResult) {
				
				if (limit != null) {
					Long previousLimit = limit.bp;

					for (RegionContent regCont : dataResult.getContents()) {

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

			queueManager.addDataRequest(
					cytobandDataSource, new DataRequest(new Region(0l, Long.MAX_VALUE, bpRegion.start.chr), 
							new HashSet<DataType>(Arrays.asList(new DataType[] {DataType.VALUE })), 
							new DataStatus()), null);
		}
	}

	public BpCoord getLimit() {
		if (limit != null && limit.bp != 0) {
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
