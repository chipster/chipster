package fi.csc.microarray.databeans;

import java.util.LinkedList;
import java.util.List;



/**
 * DataFolder is used to manage DataBean objects.
 * 
 * @see Dataset
 * @author Aleksi Kallio, hupponen
 *
 */
public class DataFolder extends DataItemBase {


	public DataFolder(DataManager manager, String name) {
		this.manager = manager;
		this.name = name;
	}
	
	private List<DataItem> children = new LinkedList<DataItem>();
	private DataManager manager;
	
	public void addChild(DataItem child) {

		// was it already connected?
		boolean wasConnected = child.getParent() != null;
	
		// connect to this
		child.setParent(this);
		
		// add
		children.add(child);
		
		// dispatch events if needed
		if (!wasConnected) {
			manager.dispatchEvent(new DataItemCreatedEvent(child));
		}
	}

	public void removeChild(DataItem child) {
		// remove connections
		child.setParent(null);
		
		// remove
		children.remove(child);		

		// dispatch events
		manager.dispatchEvent(new DataItemRemovedEvent(child));
	}

	public Iterable<DataItem> getChildren() {
		return children;
	}


	public int getChildCount() {
		return children.size();
	}

	
	public DataFolder getChildFolder(String name) {
		for (DataItem child : getChildren()) {
			if (child instanceof DataFolder && child.getName().equals(name)) {
				return (DataFolder)child;
			}
		}
		return null;
	}
	
	public String toString() {
		return getName();
	}

}
