package fi.csc.microarray.wizard.affymetrix;

import java.util.Iterator;
import java.util.LinkedList;
import javax.swing.*;

public class SelectionSet {
	LinkedList<SelectionItem> selectionList = new LinkedList<SelectionItem>();

	public void add(String publicValue, String internalValue) {
		selectionList.add(new SelectionItem(publicValue, internalValue));
	}
	
	public void initComboBox(JComboBox comboBox) {
		Iterator<SelectionItem> iter = selectionList.iterator();
		
		while (iter.hasNext()) {
			SelectionItem item = (SelectionItem)iter.next();
			comboBox.addItem(item.getPublicValue());
		}
	}
	
	public String getInternalValue(String publicValue) {
		Iterator<SelectionItem> iter = selectionList.iterator();
		
		while (iter.hasNext()) {
			SelectionItem item = (SelectionItem)iter.next();
			String value = item.getPublicValue();
			if (value.equals(publicValue))
				return item.getInternalValue();
		}
		return null;
	}
	
	public Iterator<SelectionItem> getIterator() {
		return selectionList.iterator();
	}
}
