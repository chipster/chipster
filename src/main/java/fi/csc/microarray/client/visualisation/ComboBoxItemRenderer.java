package fi.csc.microarray.client.visualisation;

import java.awt.Component;
import java.awt.Image;

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
			ImageIcon menuIcon = null;
			
			if (index < 0) {
				// small icon for closed combo box
				menuIcon = new ImageIcon(icon.getImage().getScaledInstance(16, 16, Image.SCALE_SMOOTH));
				setBorder(null);
			} else {
				// original icon size for combo box menu
				menuIcon = icon;
				setBorder(new LineBorder(this.getBackground(), 5));
			}
			
			setIcon(menuIcon);
			setText(value.toString());
			setFont(list.getFont());
		}
		
		
		return this;
	}
}
