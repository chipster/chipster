package fi.csc.microarray.client.dialog;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.HeadlessException;
import java.awt.Toolkit;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.ImportSession;
import fi.csc.microarray.client.dataimport.ImportUtils;

/**
 * Dialog for importing data from clipboard
 * 
 * @author pklemela
 *
 */
public class ClipboardImportDialog extends JDialog implements ActionListener, CaretListener{

	private final Dimension BUTTON_SIZE = new Dimension(70, 25);
	
	private static final String SKIP_TEXT = "Import directly if possible";
	private JCheckBox skipCheckBox;

	private JLabel label;
	private JButton okButton;
	private JButton cancelButton;
	private JTextField nameField;
	private JComboBox folderNameCombo;
	private ClientApplication client;

	public ClipboardImportDialog(SwingClientApplication client) {
		super(client.getMainFrame(), true);

		this.client = client;
		this.setTitle("Import data from clipboard");
		this.setModal(true);

		label = new JLabel("Filename");
		nameField = new JTextField(30);
		nameField.setText("clipboard.txt");
		nameField.addCaretListener(this);
		skipCheckBox = new JCheckBox(SKIP_TEXT);
		skipCheckBox.setSelected(true);


		folderNameCombo = new JComboBox(ImportUtils.getFolderNames(true).toArray());
		folderNameCombo.setEditable(true);
		

		okButton = new JButton("OK");
		okButton.setPreferredSize(BUTTON_SIZE);
		okButton.addActionListener(this);

		cancelButton = new JButton("Cancel");
		cancelButton.setPreferredSize(BUTTON_SIZE);
		cancelButton.addActionListener(this);			

		GridBagConstraints c = new GridBagConstraints();

		this.setLayout(new GridBagLayout());
		// Label
		c.anchor = GridBagConstraints.WEST;
		c.insets.set(10, 10, 5, 10);
		c.gridx = 0; 
		c.gridy = 0;
		this.add(label, c);

		// Combobox
		c.insets.set(0,10,10,10);		
		c.gridy++;
		this.add(nameField, c);

		c.insets.set(10,10,5,10);
		c.gridy++;
		this.add(new JLabel("Insert in folder"),c);

		c.fill = GridBagConstraints.HORIZONTAL;
		c.insets.set(0,10,10,10);
		c.gridy++;
		this.add(folderNameCombo,c);

		c.insets.set(10, 10, 10, 10);
		c.anchor = GridBagConstraints.EAST;
		c.gridy++;
		this.add(skipCheckBox, c);
		c.fill = GridBagConstraints.NONE;
		// Buttons
		c.insets.set(10, 10, 10, 10);
		c.anchor = GridBagConstraints.EAST;
		c.gridy++;
		JPanel keepButtonsRightPanel = new JPanel();
		keepButtonsRightPanel.add(okButton);
		keepButtonsRightPanel.add(cancelButton);
		this.add(keepButtonsRightPanel, c);
		
		if(this.isDataAvailable()){

			this.pack();
			this.setLocationRelativeTo(client.getMainFrame());
			this.setVisible(true);
		} else {
			JOptionPane.showMessageDialog(this, "There is no text content on the clipboard");
			this.dispose();
		}
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == okButton){
			
			File file;
			try {
				if (ImportUtils.convertToDatasetName(this.getFileName()).length() < 3){
					nameField.setText("clipboard_" + this.getFileName());
				}
				
				file = ImportUtils.createTempFile(
						ImportUtils.convertToDatasetName(this.getFileName()),
						ImportUtils.getExtension(this.getFileName()));
				
				if (pasteToFile(file, this)){	//Only if success
					ImportSession importSession = new ImportSession(ImportSession.Source.CLIPBOARD, new File[] { file }, folderNameCombo.getSelectedItem().toString(), true);
					ImportUtils.executeImport(importSession);
					this.dispose();
				}
				
			} catch (IOException ioe) {
				IOException newException = new IOException("Can't create a temporary file for the clipboard paste (reason: " + ioe.getMessage() + ").");
				client.reportException(newException);
			}			
		} else if (e.getSource() == cancelButton){
			this.dispose();
		}
	}

	public String getSelectedFolderName(){
		return folderNameCombo.getSelectedItem().toString();
	}

	public String getFileName(){
		return nameField.getText();
	}

	/**
	 * With this method, the dialog listens to interactions with the text
	 * components. The main purpose is to disable OK button if one of the
	 * inputs is empty (zero length).
	 */
	public void caretUpdate(CaretEvent e) {
		this.updateEnabledStatus();
	}

	private void updateEnabledStatus(){
		String dataSetName = nameField.getText();
		if (dataSetName.length() > 0) {
			okButton.setEnabled(true);
		} else {
			okButton.setEnabled(false);
		}
	}
	
	private boolean isDataAvailable(){
		try {
			Toolkit.getDefaultToolkit().getSystemClipboard().
			getData(DataFlavor.stringFlavor);
			return true;
			
		} catch (HeadlessException e) {
		} catch (UnsupportedFlavorException e) {
		} catch (IOException e) {
		}		
		return false;
	}

	/**
	 * Writes data from clipboard to given file
	 * 
	 * @param file file where the data is written to
	 * @return <code>true</code> if the data has been successfully written to file  
	 */
	private static boolean pasteToFile(File file, Component parentComponent){
		try{
			BufferedWriter writer = new BufferedWriter(new FileWriter(file));

			writer.append(
					(String)Toolkit.getDefaultToolkit().getSystemClipboard().
					getData(DataFlavor.stringFlavor));
			writer.close();
			return true;

		} catch(IOException ioe){
			JOptionPane.showMessageDialog(parentComponent, "Error occured while retrieving data from clipboard", "Clipboard error", JOptionPane.ERROR_MESSAGE);
			return false;
		} catch(UnsupportedFlavorException ufe){
			JOptionPane.showMessageDialog(parentComponent, "Requested data is not available", "Clipboard error", JOptionPane.ERROR_MESSAGE);
			return false;
		} catch(IllegalStateException ise){
			JOptionPane.showMessageDialog(parentComponent, "Clipboard is currently unavailable", "Clipboard error", JOptionPane.ERROR_MESSAGE);
			return false;
		}
	}
}
