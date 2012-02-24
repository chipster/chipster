package fi.csc.microarray.databeans;

import java.util.LinkedList;
import java.util.List;

/**
 * DataFolder is used to manage DataBean objects.
 * 
 * @see DataBean
 * @author Aleksi Kallio, hupponen
 * 
 */
public class DataFolder extends DataItemBase {

	public DataFolder(DataManager manager, String name) {
		this.name = name;
	}

	List<DataItem> children = new LinkedList<DataItem>(); // DataManager manages these, so package access is needed

	public Iterable<DataItem> getChildren() {
		return children;
	}

	public int getChildCount() {
		return children.size();
	}

	public DataFolder getChildFolder(String name) {
		for (DataItem child : getChildren()) {
			if (child instanceof DataFolder && child.getName().equals(name)) {
				return (DataFolder) child;
			}
		}
		return null;
	}

	public String toString() {
		return getName();
	}

}
