package fi.csc.microarray.client.dataimport.table;

import java.awt.Component;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import javax.swing.DefaultCellEditor;
import javax.swing.JTable;



/**
 * Short but essential class to enable editing of components placed in the header of 
 * JTable.
 * 
 * @author klemela
 */
public class PanelTableHeaderEditor  extends DefaultCellEditor implements ItemListener {

	PanelTableHeaderRenderer panel;

	public PanelTableHeaderEditor(PanelTableHeaderRenderer panel) {
		super(panel.getCombo());
		this.panel = panel;

		panel.getCombo().addItemListener(this);
	}

	public Component getTableCellEditorComponent(JTable table, Object value,
			boolean isSelected, int row, int column) {
		if (value instanceof Integer) {
			panel.setSelectedIndex(((Integer)value).intValue());
		} 
		return panel;
	}

	public Object getCellEditorValue() {
		return new Integer(panel.getSelectedIndex());
	}

	public void itemStateChanged(ItemEvent e) {
		super.fireEditingStopped();
	}
}

