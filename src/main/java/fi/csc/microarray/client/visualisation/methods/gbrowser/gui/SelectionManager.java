package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Selectable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;

/**
 * SelectionManager for GenomeBrowser selections. Each data has its own selections and therefore setting a
 * new selection for one of the file doesn't affect the selections of the other files. Internally the 
 * selections are handled as Selectable objects, but also row numbers (excluding comment lines) can be used 
 * to integrate selections with other visualizations.
 * 
 * @author klemela
 */
public class SelectionManager {
	
	public class SelectableSet extends HashSet<Selectable> {		
	}
	
	public class RowSet extends HashSet<Integer> {

		private Object source;

		public void setEventSource(Object source) {
			this.source = source;
		}

		public Object getEventSource() {
			return source;
		}		
	}
		
	HashMap<DataUrl, SelectableSet> selection = new HashMap<>();
	private HashMap<DataUrl, RowSet> rowSelection = new HashMap<>();
	
	private LinkedList<BrowserSelectionListener> listeners = new LinkedList<>();
	
	private GBrowser browser;
	
	public SelectionManager(GBrowser browser) {
		
		this.browser = browser;
	}

	private void add(DataUrl data, Selectable selectable, Object eventSource) {	
		
		boolean added = getSelectableSet(data).add(selectable);
		
		if (added) {
			callSelectionListeners(data, selectable, eventSource);			
		}
	}

	public SelectableSet getSelectableSet(DataUrl data) {
		if (!selection.containsKey(data)) {
			selection.put(data, new SelectableSet());
		}
		return selection.get(data);
	}
	
	private RowSet getRowSet(DataUrl data) {
		if (!rowSelection.containsKey(data)) {
			rowSelection.put(data, new RowSet());
		}
		return rowSelection.get(data);
	}

	public boolean isSelected(DataUrl data, Selectable selectable) {
		
		if (data == null) {
			return false;
		}
		
		SelectableSet set = getSelectableSet(data);
				
		boolean internal =  set.contains(selectable);		
		boolean external = false;
		
		if (selectable.getIndexKey() != null && selectable.getIndexKey().getRowNumber() != null) {
			
			Integer row = (Integer)(int)selectable.getIndexKey().getRowNumber();
			
			if (row != null && rowSelection.containsKey(data)) {
				external = rowSelection.get(data).contains(row);
			} else {
				external = false;
			}
			
			if (external && !internal) {
				//Source information tells if this the selection change originated from this browser 
				//or other visualizations. This information is needed later to decide whether the change
				//should be reported to other visualizations or not.
				Object eventSource = getRowSet(data).getEventSource();
				add(data, selectable, eventSource);
				getRowSet(data).remove(row);
			}
		}
		
		return internal || external;
	}
	
	public void set(DataUrl data, Selectable selectable) {
		
		getSelectableSet(data).clear();
		getRowSet(data).clear();
		
		if (selectable != null) {
			this.add(data, selectable, browser);
		} else {
			//Report that selection was cleared
			callSelectionListeners(data, null, browser);
		}
	}

	private void callSelectionListeners(DataUrl data, Selectable selectable, Object source) {
		for (BrowserSelectionListener listener : listeners) {
			listener.selectionChanged(data, selectable, source);
		}
	}

	public void addSelectionListener(BrowserSelectionListener listener) {
		listeners.add(listener);
	}

	public String getSelectionText() {
		StringBuilder builder = new StringBuilder();
		
		LinkedList<DataUrl> datasIntrackOrder = new LinkedList<>();
		
		//Use track order to print selections
		for (Track track : browser.getPlot().getDataView().getTracks()) {
			DataUrl data = track.getDataUrl();

			if (data != null) {
				//Many tracks may use the same data, but it is enough to print the selection once
				if (!datasIntrackOrder.contains(data)) {
					datasIntrackOrder.add(data);
				}
			}
		}
		
		for (DataUrl data : datasIntrackOrder) {
			SelectableSet set = selection.get(data);
			if (set != null) {
				for (Selectable item : selection.get(data)) {					
					builder.append(item.getText());
					builder.append("\n\n");
				}
			}
		}
		return builder.toString();
	}

	public void toggle(DataUrl data, Selectable selectable) {
		SelectableSet set = getSelectableSet(data);
		if (set.contains(selectable)) {
			remove(data, selectable, browser);
		} else {
			add(data, selectable, browser);
		}
	}

	private void remove(DataUrl data, Selectable selectable, Object eventSource) {
		boolean removed = getSelectableSet(data).remove(selectable);
		if (removed) {
			callSelectionListeners(data, selectable, eventSource);
		}
	}
	
	/**
	 * @return all dataUrls that have at least one Selectable selected at the moment. 
	 */
	public Set<DataUrl> getDataUrls() {
		return selection.keySet();		
	}
	
	/**
	 * SelectionManager is based on Selectable objects, which we don't have. Just clear any existing selections,
	 * store the row numbers and redraw. Later in method isSelected the Selectable objects are compared with
	 * stored row numbers.
	 * 
	 * @param data
	 * @param eventSource
	 */
	public void setRowSelections(DataUrl data, HashSet<Integer> rows, Object eventSource) {
		
		//Clear all existing selections for this data
		getSelectableSet(data).clear();
		
		RowSet rowSet = getRowSet(data);		
		rowSet.clear();
		rowSet.setEventSource(eventSource);		
		rowSet.addAll(rows);
		
		if (SelectionManager.this.browser.getPlot() != null) {
			SelectionManager.this.browser.getPlot().redraw();
		}
	}
}
