package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;
import java.io.File;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;


public class FsfStatus {
	
	public boolean poison; // All threads should send this forward and end themselves
	public long areaRequestCount;
	public long fileRequestCount;
	public long fileResultCount;
//	public long fileThreadWait;
//	public long treeThreadWait;
	public boolean clearQueues;
	public boolean concise;
	public boolean debug;
	
	private Set<Collection> clearedAlready = new HashSet<Collection>();
	public File file;
	
	public void maybeClearQueue(Collection queue){
		if(clearQueues && !clearedAlready.contains(queue)){

			clearedAlready.add(queue);
			queue.clear();			
		}
	}
}
