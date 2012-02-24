package fi.csc.microarray.databeans;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.databeans.features.BoolFalseFeature;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.databeans.features.table.TableColumnProvider.MatrixParseSettings;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.IOUtils;

/**
 * DataFolder is used to manage DataBean objects.
 * 
 * @see DataBean
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

		if (child instanceof DataBean) {
			try {
				manager.doBackwardsCompatibleTypeTagInitialisation((DataBean) child);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

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
				return (DataFolder) child;
			}
		}
		return null;
	}

	public String toString() {
		return getName();
	}

}
