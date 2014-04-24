package fi.csc.microarray.client.visualisation;

import java.awt.Component;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.ListCellRenderer;
import javax.swing.border.LineBorder;

class ComboBoxRenderer extends JLabel implements ListCellRenderer<VisualisationMethod> {
	public ComboBoxRenderer() {
		setOpaque(true);
		setHorizontalAlignment(LEFT);
		setVerticalAlignment(CENTER);		
	}

	/*
	 * This method finds the image and text corresponding
	 * to the selected value and returns the label, set up
	 * to display the text and image.
	 */
	@Override
	public Component getListCellRendererComponent(
			JList<? extends VisualisationMethod> list,
			VisualisationMethod value,
			int index,
			boolean isSelected,
			boolean cellHasFocus) {
//		Get the selected index. (The index param isn't
//		always valid, so just use the value.)
		//int selectedIndex = ((Integer)value).intValue();

		if (isSelected) {
			setBackground(list.getSelectionBackground());
			setForeground(list.getSelectionForeground());
		} else {
			setBackground(list.getBackground());
			setForeground(list.getForeground());
		}

		ImageIcon icon = null;
		if (value instanceof VisualisationMethod) {
			icon = ((VisualisationMethod)value).getIcon(); 
		}
		
		if (value != null) {
						
			if (index < 0) {
				// closed menu
				setIcon(null);
				setBorder(null);
				//setText(value.toString());
			} else {
				// open menu
				setIcon(icon);
				setBorder(new LineBorder(this.getBackground(), 5));
			}
			
			// both closed and open menu look better when there is some empty space in front of the text
			setText("  " + value.toString());
			setFont(list.getFont());
		}
		
		
		return this;
	}
}
