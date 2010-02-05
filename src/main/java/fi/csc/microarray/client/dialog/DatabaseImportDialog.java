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

import javax.swing.BorderFactory;
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
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.parameter.ImportParameterPanel;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

/**
 * Dialog for importing data from clipboard
 * 
 * @author pklemela
 * 
 */
public class DatabaseImportDialog extends JDialog implements ActionListener,
		CaretListener {

	private final Dimension BUTTON_SIZE = new Dimension(70, 25);

	private JLabel titleLabel;
	private JLabel descriptionLabel;
	private JButton okButton;
	private JButton cancelButton;
	private JComboBox folderNameCombo;
	private ClientApplication application;
	private Operation operation;

	public DatabaseImportDialog(SwingClientApplication application, String databaseName, Operation operation) throws MicroarrayException {
		super(application.getMainFrame(), true);

		this.application = application;
		this.operation = operation;
		this.setTitle("Import from " + databaseName);
		this.setModal(true);
		//this.setPreferredSize(new Dimension(500, 400));

		titleLabel = new JLabel("<html><p style=\"font-weight:bold;font-size:120%\">Import data from " + databaseName + "</p></html>");
		//titleLabel.setPreferredSize(new Dimension(300, 300));
		descriptionLabel = new JLabel("<html>" + operation.getDescription() + "</html>");
		

		folderNameCombo = new JComboBox(ImportUtils.getFolderNames(true)
				.toArray());
		folderNameCombo.setEditable(true);

		okButton = new JButton("Import");
		okButton.setPreferredSize(BUTTON_SIZE);
		okButton.addActionListener(this);

		cancelButton = new JButton("Cancel");
		cancelButton.setPreferredSize(BUTTON_SIZE);
		cancelButton.addActionListener(this);

		GridBagConstraints c = new GridBagConstraints();

		this.setLayout(new GridBagLayout());
		c.weightx = 1.0;
		c.weighty = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;
//		c.gridwidth = GridBagConstraints.REMAINDER;
		
		// title label
		c.anchor = GridBagConstraints.NORTHWEST;
		c.insets.set(10, 10, 5, 10);
		c.gridx = 0;
		c.gridy = 0;
//		titleLabel.setBorder(BorderFactory.createMatteBorder(1,1,1,1,VisualConstants.OPERATION_LIST_BORDER_COLOR));
		this.add(titleLabel, c);

		// description label
		c.anchor = GridBagConstraints.NORTHWEST;
		c.gridy++;
//		descriptionLabel.setBorder(BorderFactory.createMatteBorder(1,1,1,1,VisualConstants.OPERATION_LIST_BORDER_COLOR));
		this.add(descriptionLabel, c);
		
		
		// parameter panel
		c.weighty = 120.0;
		c.fill = GridBagConstraints.BOTH;
		c.insets.set(0, 10, 10, 10);
		c.gridy++;
		ImportParameterPanel parameterPanel = new ImportParameterPanel(operation, null);
		this.add(parameterPanel, c);

		c.fill = GridBagConstraints.NONE;
		c.insets.set(10, 10, 5, 10);
		c.gridy++;
//		this.add(new JLabel("Insert in folder"), c);

		c.fill = GridBagConstraints.NONE;
		c.insets.set(0, 10, 10, 10);
		c.gridx = GridBagConstraints.RELATIVE;
//		this.add(folderNameCombo, c);

		JPanel folderPanel = new JPanel();
		folderPanel.add(new JLabel("Insert in folder"));
		folderPanel.add(folderNameCombo);
		
		this.add(folderPanel, c);
		
		c.insets.set(10, 10, 10, 10);
		c.anchor = GridBagConstraints.SOUTHEAST;
		c.gridx = 0;
		c.gridy++;
		c.fill = GridBagConstraints.NONE;

		
		// Buttons
		c.insets.set(10, 10, 10, 10);
		c.anchor = GridBagConstraints.SOUTHEAST;
		c.gridy++;
		c.weightx = 1.0;
		c.weighty = 20.0;
		c.fill = GridBagConstraints.NONE;
		
		JPanel keepButtonsRightPanel = new JPanel(new GridBagLayout());
		keepButtonsRightPanel.setBorder(BorderFactory.createMatteBorder(1,1,1,1,VisualConstants.OPERATION_LIST_BORDER_COLOR));
		GridBagConstraints buttonConstraints = new GridBagConstraints();
		
		buttonConstraints.weightx = 1.0;
		buttonConstraints.weighty = 1.0;
		buttonConstraints.anchor = GridBagConstraints.EAST;
		buttonConstraints.insets.set(10, 10, 10, 10);
		keepButtonsRightPanel.add(cancelButton, buttonConstraints);
		buttonConstraints.gridx = GridBagConstraints.RELATIVE;
		keepButtonsRightPanel.add(okButton, buttonConstraints);
		this.add(keepButtonsRightPanel, c);

		if (this.isDataAvailable()) {

			this.pack();
			this.setLocationRelativeTo(application.getMainFrame());
			this.setVisible(true);
		} else {
			JOptionPane.showMessageDialog(this,
					"There is no text content on the clipboard");
			this.dispose();
		}
	}

	public void actionPerformed(ActionEvent e) {

		// start the import task
		if (e.getSource() == okButton) {
			try {
				application.executeOperation(operation);
			} catch (Exception me) {
				application.reportException(me);
			} finally {
				this.dispose();
			}
		} 
		
		// cancel import
		else if (e.getSource() == cancelButton) {
			this.dispose();
		}
	}

	public String getSelectedFolderName() {
		return folderNameCombo.getSelectedItem().toString();
	}

	/**
	 * With this method, the dialog listens to interactions with the text
	 * components. The main purpose is to disable OK button if one of the inputs
	 * is empty (zero length).
	 */
	public void caretUpdate(CaretEvent e) {
		this.updateEnabledStatus();
	}

	private void updateEnabledStatus() {
//		String dataSetName = nameField.getText();
//		if (dataSetName.length() > 0) {
//			okButton.setEnabled(true);
//		} else {
//			okButton.setEnabled(false);
//		}
	}

	private boolean isDataAvailable() {
		try {
			Toolkit.getDefaultToolkit().getSystemClipboard().getData(
					DataFlavor.stringFlavor);
			return true;

		} catch (HeadlessException e) {
		} catch (UnsupportedFlavorException e) {
		} catch (IOException e) {
		}
		return false;
	}
}
