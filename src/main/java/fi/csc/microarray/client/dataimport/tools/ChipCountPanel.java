package fi.csc.microarray.client.dataimport.tools;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JCheckBox;
import javax.swing.JLabel;

import org.jdesktop.swingx.JXTaskPane;

import fi.csc.microarray.client.dataimport.ColumnType;
import fi.csc.microarray.client.dataimport.ImportScreen;

public class ChipCountPanel extends JXTaskPane implements ActionListener{
	
	private JCheckBox[] keyColumnCheckBoxes;
	private JLabel[] keyColumnCounterLabels;
	private JLabel[] keyColumnTotalCountLabels;
	private JLabel firstRowLabel;
	
	private ImportScreen screen;
	
	/**
	 * Key values to keep them in nice order
	 */
	protected static final ColumnType[] KEYS =
		new ColumnType[] {
			ColumnType.SAMPLE_LABEL,
			ColumnType.SAMPLE_BG_LABEL,
			ColumnType.CONTROL_LABEL,
			ColumnType.CONTROL_BG_LABEL,
			ColumnType.FLAG_LABEL,
			ColumnType.IDENTIFIER_LABEL,
			ColumnType.ANNOTATION_LABEL
			};
	
	
	public static boolean isCountedKey(ColumnType key) {
		for (int i = 0; i < KEYS.length; i++) {
			if (KEYS[i].equals(key)) {
				return true;
			}
		}
		return false;
	}
	
	public ChipCountPanel(ImportScreen screen) {
		this.setLayout(new GridBagLayout());		
		this.screen = screen;
				
		keyColumnCheckBoxes = new JCheckBox[KEYS.length];
		keyColumnCounterLabels = new JLabel[KEYS.length];
		keyColumnTotalCountLabels = new JLabel[KEYS.length];
		
		for (int i = 0; i < KEYS.length; i++) {
			if (i == 0) {
				// Name label added separately to keep it black
				keyColumnCheckBoxes[i] = new JCheckBox();
				// To place the label same way as checkboxes naturally do
				keyColumnCheckBoxes[i].setPreferredSize(new Dimension(17, 19));
				keyColumnCheckBoxes[i].setSelected(true);  // always at least one sample column..
				keyColumnCheckBoxes[i].setEnabled(false);  // ..so this can't be unchecked
			} else {
				keyColumnCheckBoxes[i] = new JCheckBox(KEYS[i].toString());
				keyColumnCheckBoxes[i].setSelected(true);
				keyColumnCheckBoxes[i].setActionCommand(KEYS[i].toString());
				keyColumnCheckBoxes[i].addActionListener(this);
			}
			
			keyColumnCounterLabels[i] = new JLabel("0");
			keyColumnCounterLabels[i].setOpaque(true);
			
			keyColumnTotalCountLabels[i] = new JLabel(" of 1");
			keyColumnTotalCountLabels[i].setOpaque(true);
			
			Color rowBg = Color.white;//ImportPreviewTable.getColorForLabel(KEY_NAMES[i]);
			keyColumnCheckBoxes[i].setBackground(rowBg);
			keyColumnCounterLabels[i].setBackground(rowBg);
			keyColumnTotalCountLabels[i].setBackground(rowBg);
			
			/*
			if (i > 0) {
				Border topLine = VisualConstants.BORDER_DARKLINE_TOP;
				keyColumnCheckBoxes[i].setBorder(topLine);
				keyColumnCounters[i].setBorder(topLine);
				keyColumnTotalCounters[i].setBorder(topLine);
			}
			*/
		}
		
		firstRowLabel = new JLabel(KEYS[0].toString());
		firstRowLabel.setOpaque(true);
		firstRowLabel.setBackground(Color.white);//ImportPreviewTable.getColorForLabel(KEY_NAMES[0]));
		
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0; 
		c.gridy = 0;
		c.fill = GridBagConstraints.BOTH;
		
		for (int i = 0; i < KEYS.length; i++) {
			c.gridx = 0; c.gridy++;
			c.anchor = GridBagConstraints.WEST;			
			
			if (i == 0) {
				//To have black text with disabled checkbox				
				c.weightx = 0;
				c.gridwidth = 1;
				this.add(keyColumnCheckBoxes[0], c);
				c.gridx = 1;
				c.weightx = 1;				
				c.gridwidth = 2;
				this.add(firstRowLabel, c);
			} else {
				c.gridwidth = 3;
				this.add(keyColumnCheckBoxes[i], c);
			}
			c.gridx = 3;
			c.weightx = 0;
			c.gridwidth = 1;
			c.anchor = GridBagConstraints.EAST;
			this.add(keyColumnCounterLabels[i], c);
			c.gridx = 4;
			this.add(keyColumnTotalCountLabels[i], c);
		}

		this.setTitle("Channels");
		this.setExpanded(false);
	}
	
	public void updateKeyColumnCounter(int keyNumber) {
		
		int countOfProperlySetKeys = screen.getColumnTypeManager().getCountOfCorrectlySet(KEYS[keyNumber]);
		int countOfKeys = screen.getColumnTypeManager().getCountOfType(KEYS[keyNumber]);
		
		keyColumnCounterLabels[keyNumber].setText(""+countOfKeys);
		
		if (keyColumnCheckBoxes[keyNumber].isSelected()) {
			if(KEYS[keyNumber].equals(ColumnType.IDENTIFIER_LABEL)){
				keyColumnTotalCountLabels[keyNumber].setText(" of 1");
			} else if(KEYS[keyNumber].equals(ColumnType.ANNOTATION_LABEL)){
				keyColumnTotalCountLabels[keyNumber].setText("");
			} else {
				keyColumnTotalCountLabels[keyNumber].setText(
						" of " + screen.getColumnTypeManager().getCountOfType(ColumnType.SAMPLE_LABEL));
			}
		} else {
			keyColumnTotalCountLabels[keyNumber].setText(
					" of " + 0);
		}

		int sampleCount = screen.getColumnTypeManager().getCountOfType(ColumnType.SAMPLE_LABEL);
		
		// Color the label red if more keys selected than chips exists
		if(KEYS[keyNumber].equals(ColumnType.ANNOTATION_LABEL)){
			keyColumnCounterLabels[keyNumber].setForeground(Color.BLACK);
			keyColumnTotalCountLabels[keyNumber].setForeground(Color.BLACK);
		} else if(KEYS[keyNumber].equals(ColumnType.IDENTIFIER_LABEL)){
			if(countOfKeys > 1 || countOfKeys > sampleCount){
				keyColumnCounterLabels[keyNumber].setForeground(Color.RED);
				keyColumnTotalCountLabels[keyNumber].setForeground(Color.RED);
			} else {
				keyColumnCounterLabels[keyNumber].setForeground(Color.BLACK);
				keyColumnTotalCountLabels[keyNumber].setForeground(Color.BLACK);
			}
		} else {
			if(countOfProperlySetKeys == sampleCount 
					&& countOfKeys == sampleCount && sampleCount > 0){
				keyColumnCounterLabels[keyNumber].setForeground(Color.GREEN.darker());
				keyColumnTotalCountLabels[keyNumber].setForeground(Color.GREEN.darker());
			} else if(countOfKeys > sampleCount || (countOfKeys == sampleCount && countOfProperlySetKeys < sampleCount)){
				keyColumnCounterLabels[keyNumber].setForeground(Color.RED);
				keyColumnTotalCountLabels[keyNumber].setForeground(Color.RED);
			} else {
				keyColumnCounterLabels[keyNumber].setForeground(Color.BLACK);
				keyColumnTotalCountLabels[keyNumber].setForeground(Color.BLACK);
			}
		}
		
		this.repaint();
	}
	
	public void updateAllKeyColumnCounters() {
		for (int i = 0; i < KEYS.length; i++) {
			this.updateKeyColumnCounter(i);
		}
	}
	
	public void setEnabled(boolean enabled) {
		for (int i = 0; i < KEYS.length; i++) {
			firstRowLabel.setEnabled(enabled);
			if (i > 0) {
				keyColumnCheckBoxes[i].setEnabled(enabled);
			}
			keyColumnCounterLabels[i].setEnabled(enabled);
			keyColumnTotalCountLabels[i].setEnabled(enabled);
		}
	}
	
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		if (source instanceof JCheckBox) {
			for (int i = 0; i < keyColumnCheckBoxes.length; i++) {
				if (keyColumnCheckBoxes[i] == source) {
					this.updateKeyColumnCounter(i);
					boolean isSelected = keyColumnCheckBoxes[i].isSelected();
					keyColumnCounterLabels[i].setEnabled(isSelected);
					keyColumnTotalCountLabels[i].setEnabled(isSelected);
					break;
				}
			}
		}
	}
}
