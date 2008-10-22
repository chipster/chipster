package fi.csc.microarray.client.dataimport.tools;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.jdesktop.swingx.JXTaskPane;

import fi.csc.microarray.client.dataimport.ImportScreen;
import fi.csc.microarray.client.dataimport.events.ColumnTitlesChangedEvent;
import fi.csc.microarray.client.dataimport.events.ConversionModelChangeListener;
import fi.csc.microarray.client.dataimport.events.DecimalSeparatorChangedEvent;
import fi.csc.microarray.client.dataimport.events.DelimiterChangedEvent;
import fi.csc.microarray.client.dataimport.events.FooterChangedEvent;
import fi.csc.microarray.client.dataimport.events.HeaderChangedEvent;
import fi.csc.microarray.client.dataimport.events.InputFileChangedEvent;
import fi.csc.microarray.client.dataimport.events.TitleRowChangedEvent;
import fi.csc.microarray.client.dataimport.table.ImportPreviewTable;
import fi.csc.microarray.client.dataimport.trimmer.DataTrimmingOperation;
import fi.csc.microarray.client.dataimport.trimmer.NormalStringReplace;
import fi.csc.microarray.client.dataimport.trimmer.ReqularExpressionStringReplace;

/**
 * Find-and-replace functionality for the import screen.
 */
public class DataTrimmingPanel extends JXTaskPane implements ActionListener, ListSelectionListener, ConversionModelChangeListener {

	private static final String ALL_COLUMNS = "All";
	
	private JLabel columnLabel;
	private JLabel lookForLabel;
	private JLabel replaceWithLabel;
	private JTextField lookForField;
	private JTextField replaceWithField;
	private JCheckBox useRegexpsCheckBox;
	private JButton replaceButton;
	private JButton undoReplaceButton;
	private JComboBox columnCombobox;
	
	private ImportScreen screen;
	private ImportPreviewTable table;
	
	public DataTrimmingPanel(ImportScreen screen) {
		this.setLayout(new GridBagLayout());
		
		this.screen = screen;
		this.table = screen.getTableFrame().getTable();
		
		columnCombobox = new JComboBox();
		this.updateComboboxItems();
		
		columnLabel = new JLabel("Column:");
		lookForLabel = new JLabel("Look For:");
		replaceWithLabel = new JLabel("Replace With:");
		
		lookForField = new JTextField();
		lookForField.addCaretListener(new CaretListener() {
			public void caretUpdate(CaretEvent e) {
				if (lookForField.getText().length() == 0) {
					replaceButton.setEnabled(false);
				} else {
					replaceButton.setEnabled(true);
				}
			}
		});
		
		replaceWithField = new JTextField();
		
		useRegexpsCheckBox = new JCheckBox("Use Regular Expressions");
		useRegexpsCheckBox.setPreferredSize(new Dimension(165, 18));
		useRegexpsCheckBox.setSelected(false);
		useRegexpsCheckBox.setOpaque(false);
		useRegexpsCheckBox.addActionListener(this);
		
		replaceButton = new JButton("Replace");
		replaceButton.addActionListener(this);
		replaceButton.setEnabled(false);
		
		undoReplaceButton = new JButton("Undo");
		undoReplaceButton.addActionListener(this);
		undoReplaceButton.setEnabled(false);
		
		GridBagConstraints c = new GridBagConstraints();
		
		c = new GridBagConstraints();
		c.gridx = 0; c.gridy = 0;
		c.weightx = 1;
		c.gridwidth = 2;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;

		this.add(columnLabel, c);
		
		c.gridy++;
		this.add(columnCombobox, c);		
		c.gridy++;
		this.add(lookForLabel, c);		
		c.gridy++;
		this.add(lookForField, c);		
		c.gridy++;
		this.add(replaceWithLabel, c);		
		c.gridy++;
		this.add(replaceWithField, c);		
		c.gridy++;
		this.add(useRegexpsCheckBox, c);		
		c.gridwidth = 1;
		c.gridy++;
		this.add(replaceButton, c);		
		c.gridx++;
		this.add(undoReplaceButton, c);

		this.setTitle("Data Modification");
		//this.setExpanded(false);
	}
	
	private void updateComboboxItems() {
		String[] columnTitles = screen.getConversionModel().getColumnTitles();
		
		// The length should be columnTitles.length -1 because the first column 
		// is ignored, but because the "All" option is added, the lenght is same 
		// as columnTitles length
		String[] columnTitlesNotTheFirst = new String[columnTitles.length];
		for(int i = 0; i < columnTitlesNotTheFirst.length-1; i++){
			columnTitlesNotTheFirst[i] = columnTitles[i+1];
		}
		columnTitlesNotTheFirst[columnTitles.length-1] = ALL_COLUMNS;
		
		columnCombobox.setModel(new DefaultComboBoxModel(columnTitlesNotTheFirst));
	}

	public void actionPerformed(ActionEvent e) {
		
		// Replace
		if(e.getSource() == replaceButton){
			
			boolean allSelected = false;
			int columnIndex = -1;
			if(this.columnCombobox.getSelectedItem().toString().equals(ALL_COLUMNS)){
				allSelected = true;
			} else {
				columnIndex = screen.getTableFrame().getTable().getColumnModel().getColumnIndex(this.columnCombobox.getSelectedItem().toString());
				allSelected = false;
			}
			
			if(useRegexpsCheckBox.isSelected()){	
				
				if(allSelected){
					// Ignore the first
					screen.getDataTrimmer().pushOperation(new ReqularExpressionStringReplace(lookForField.getText(), replaceWithField.getText(), DataTrimmingOperation.ALL_COLUMNS));
				} else {
					screen.getDataTrimmer().pushOperation(new ReqularExpressionStringReplace(lookForField.getText(), replaceWithField.getText(), columnIndex));
				}
				
				table.repaint();
			} else {
				if(allSelected){
					screen.getDataTrimmer().pushOperation(new NormalStringReplace(lookForField.getText(), replaceWithField.getText(), false, DataTrimmingOperation.ALL_COLUMNS));
				} else {
					screen.getDataTrimmer().pushOperation(new NormalStringReplace(lookForField.getText(), replaceWithField.getText(), false, columnIndex));
				}
			}
			undoReplaceButton.setEnabled(true);
			table.repaint();
		} 
		
		// Undo
		else if(e.getSource() == undoReplaceButton){
			screen.getDataTrimmer().popOperation();
			if(screen.getDataTrimmer().getOperationCount() == 0){
				undoReplaceButton.setEnabled(false);
			} else {
				undoReplaceButton.setEnabled(true);
			}
			table.repaint();
		}
	}

	public void valueChanged(ListSelectionEvent e) {
		this.setEnabled(true);
	}

	public void decimalSeparatorChanged(DecimalSeparatorChangedEvent e) {
		// Do nothign
		
	}

	public void delimiterChanged(DelimiterChangedEvent e) {
		// Do nothing
		
	}

	public void footerChanged(FooterChangedEvent e) {
		// Do nothing
		
	}

	public void headerChanged(HeaderChangedEvent e) {
		// Do nothing
	}

	public void titleRowChanged(TitleRowChangedEvent e) {
		updateComboboxItems();
	}

	public void columnTitlesChanged(ColumnTitlesChangedEvent e) {
		updateComboboxItems();
	}

	public void inputFileChanged(InputFileChangedEvent e) {
		// TODO Auto-generated method stub
		
	}
}
