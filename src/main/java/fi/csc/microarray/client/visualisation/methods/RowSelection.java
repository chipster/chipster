//package fi.csc.microarray.client.visualisation.methods;
//
//import java.util.HashMap;
//import java.util.Map.Entry;
//
//import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
//
//public class RowSelection {
//
//	private HashMap<DataUrl, DataRowSelection> selectionMap = new HashMap<>();
//
//	public RowSelection() {
//	}
//	
//	public void add(DataUrl data, Integer i, Object source) {
//		if (!selectionMap.containsKey(data)) {
//			selectionMap.put(data, new DataRowSelection());
//		}
//		DataRowSelection dataRowSelection = selectionMap.get(data);
//		dataRowSelection.add((Integer)i);
//		dataRowSelection.setSource(source);
//	}
//	
//	public boolean contains(DataUrl dataUrl, Integer row) {
//		DataRowSelection dataRowSelection = selectionMap.get(dataUrl);
//		if (dataRowSelection != null) {
//			return dataRowSelection.contains(row);
//		}
//		return false;
//	}
//
//	public Object getSource(DataUrl dataUrl) {
//		DataRowSelection dataRowSelection = selectionMap.get(dataUrl);
//		if (dataRowSelection != null) {
//			return dataRowSelection.getSource();
//		}
//		return null;
//	}
//
//	public void addAll(RowSelection rowSelection) {
//		for (Entry<DataUrl, DataRowSelection> entry : rowSelection.selectionMap.entrySet()) {
//			for (Integer i : entry.getValue()) {
//				this.add(entry.getKey(), i, entry.getValue().getSource());
//			}
//		}
//	}
//
//	public void clear(DataUrl dataUrl) {
//		DataRowSelection dataRowSelection = selectionMap.get(dataUrl);
//		if (dataRowSelection != null) {
//			dataRowSelection.clear();
//		}
//	}
//
//	public void remove(DataUrl dataUrl, Integer row) {
//		DataRowSelection dataRowSelection = selectionMap.get(dataUrl);
//		if (dataRowSelection != null) {
//			dataRowSelection.remove((Integer)row);
//		}
//	}
//}
