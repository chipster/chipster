package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.HashSet;
import java.util.Queue;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TreeNode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;

/**
 * <p>Fast seek file status.</p>
 * 
 * <p>A generic status field that travels with all requests. It is used to synchronise the different parts of the system.</p> 
 */
public class FsfStatus implements Cloneable {

	/**
	 * All threads should send this forward and end themselves
	 */
	public boolean poison;

	public long areaRequestCount;
	public long fileRequestCount;
	public long fileResultCount;
	public boolean clearQueues;
	public boolean concise;
	public boolean debug;

	/**
	 * All objects originating from a single area request share the same instance of this Set.
	 */
	private Set<Object> clearedAlready = new HashSet<Object>();
	
	public DataSource file;

	public TreeNode bpSearchSource;


	public void maybeClearQueue(Object requestQueue) {
		
		if (clearQueues && !clearedAlready.contains(requestQueue)) {
		
			clearedAlready.add(requestQueue);
			
			((Queue<?>) requestQueue).clear();
		}
	}
	
	@Override
	public FsfStatus clone() throws CloneNotSupportedException {
		FsfStatus status = (FsfStatus)super.clone();
		// do not clone clearedAlready, it must be shared between clones
		return status;
	}
	
}
