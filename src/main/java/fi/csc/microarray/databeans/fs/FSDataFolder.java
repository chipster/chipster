package fi.csc.microarray.databeans.fs;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataFolderBase;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataItemCreatedEvent;
import fi.csc.microarray.databeans.DataItemRemovedEvent;


public class FSDataFolder extends DataFolderBase {

	private String name;
	private DataFolder parent = null;
	private List<DataItem> children = new LinkedList<DataItem>();
	private FSDataManager manager;
	
	public FSDataFolder(FSDataManager manager, String name) {
		this.manager = manager;
		this.name = name;
	}
	
	/**
	 * @return Parent folder of this folder.
	 */
	public DataFolder getParent() {
		return parent;
	}

	private void setParent(FSDataFolder dataFolder) {
		this.parent = dataFolder;		
	}


	public void addChild(DataItem child) {

		// was it already connected?
		boolean wasConnected = child.getParent() != null;
		
		// connect to this
		if (child instanceof FSDataFolder) {
			((FSDataFolder)child).setParent(this);
		} else if (child instanceof FSDataBean) {
			((FSDataBean)child).setParent(this);
		} else {
			throw new IllegalArgumentException("cannot add child of type " + child.getClass().getSimpleName());
		}
		
		// add
		children.add(child);
		
		// dispatch events if needed
		if (!wasConnected) {
			manager.dispatchEvent(new DataItemCreatedEvent(child));
		}
	}

	public void removeChild(DataItem child) {
		// remove connections
		if (child instanceof FSDataFolder) {
			((FSDataFolder)child).setParent(null);
		} else if (child instanceof FSDataBean) {
			((FSDataBean)child).setParent(null);
		} else {
			throw new IllegalArgumentException("cannot remove child of type " + child.getClass().getSimpleName());
		}
		
		// remove
		children.remove(child);		

		// dispatch events
		manager.dispatchEvent(new DataItemRemovedEvent(child));
	}

	public Iterable<DataItem> getChildren() {
		return children;
	}

	public String getName() {
		return name; 
	}
    
    public void setName(String newName) {
        name = newName;
    }

	public int getChildCount() {
		return children.size();
	}

}
