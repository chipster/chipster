package fi.csc.microarray.wizard;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import fi.csc.microarray.client.operation.Operation.ResultListener;
import fi.csc.microarray.databeans.DataBean;

public class ResultBlocker implements ResultListener {

	private CountDownLatch latch = new CountDownLatch(1);
	private Iterable<DataBean> results;
	private int enforcedResultCount;

	public ResultBlocker() {
		this(-1);
	}

	public ResultBlocker(int enforcedResultCount) {
		this.enforcedResultCount = enforcedResultCount;
	}
	
	public void resultData(Iterable<DataBean> results) {
		this.results = results;
		latch.countDown();		
	}
	
	public void noResults() {
		this.results = null;
		latch.countDown();
	}
	
	public void blockUntilDone() {
		try {
			latch.await();
		} catch (InterruptedException e) {
			// should not happen
			throw new RuntimeException(e);
		}
	}
	
	public List<DataBean> getResults() {
		if (results == null) {
			throw new IllegalStateException("no results");
		}
		
		Iterator<DataBean> iter = results.iterator();
		LinkedList<DataBean> beans = new LinkedList<DataBean>();

		while (iter.hasNext()) {
			beans.add(iter.next());
		}
		
		if (enforcedResultCount != -1 && beans.size() != enforcedResultCount) {
			throw new IllegalStateException("there are more or less than " + enforcedResultCount + " results");
		}
		
		return beans;
	}

}
