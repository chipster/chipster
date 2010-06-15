package fi.csc.microarray.client.tasks;

import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import fi.csc.microarray.client.operation.Operation.ResultListener;
import fi.csc.microarray.databeans.Dataset;

public class ResultBlocker implements ResultListener {

	private CountDownLatch latch = new CountDownLatch(1);
	private Iterable<Dataset> results;
	private int enforcedResultCount;

	public ResultBlocker() {
		this(-1);
	}

	public ResultBlocker(int enforcedResultCount) {
		this.enforcedResultCount = enforcedResultCount;
	}
	
	public void resultData(Iterable<Dataset> results) {
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
	
	public List<Dataset> getResults() {
		if (results == null) {
			throw new IllegalStateException("no results");
		}
		
		Iterator<Dataset> iter = results.iterator();
		LinkedList<Dataset> beans = new LinkedList<Dataset>();

		while (iter.hasNext()) {
			beans.add(iter.next());
		}
		
		if (enforcedResultCount != -1 && beans.size() != enforcedResultCount) {
			throw new IllegalStateException("there are more or less than " + enforcedResultCount + " results");
		}
		
		return beans;
	}

}
