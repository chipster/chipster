package fi.csc.microarray.client.dialog;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataItem;

/**
 * A modal dialog which allows the user to rename the selected data item.
 * 
 * @author Janne Käki
 *
 */
public class RenameDialog extends JDialog implements ActionListener, CaretListener {

	private final Dimension BUTTON_SIZE = new Dimension(75, 25);
	
	private JTextField dataItemNameField;
	private JButton okButton;
	private JButton cancelButton;
	
	private ClientApplication client;
	private DataItem data;
	
	/**
	 * Creates a new data item rename dialog and shows it.
	 * 
	 * @param owner The owning frame component of the dialog to be created.
	 * 				In this case, it would be the GUI mainframe class.
	 * @param data  The dataset to be eventually renamed.
	 */
	public RenameDialog(SwingClientApplication client, DataItem data) {
		super(client.getMainFrame(),                
                data instanceof DataBean ? "Rename Data Item" : "Rename Data Folder",
                true);
		
		this.client = client;
		this.data = data;
		
		dataItemNameField = new JTextField(15);
		dataItemNameField.setText(data.getName());
        dataItemNameField.addCaretListener(this);
        dataItemNameField.addActionListener(this);
		
		okButton = new JButton("OK");
		okButton.setPreferredSize(BUTTON_SIZE);
		cancelButton = new JButton("Cancel");
		cancelButton.setPreferredSize(BUTTON_SIZE);
		
		JPanel contentPane = new JPanel(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.WEST;
		c.insets = new Insets(10, 20, 5, 10);
		c.gridx = 0; c.gridy = 0;
		c.gridwidth = 2;
		JLabel dataSetNameLabel = new JLabel("New Name for " + data.getName());
		contentPane.add(dataSetNameLabel, c);
		c.gridy = 1;
		contentPane.add(dataItemNameField, c);
		c.gridy = 2;
		c.gridwidth = 1;
		contentPane.add(okButton, c);
		c.gridx = 1;
		contentPane.add(cancelButton, c);
		this.setContentPane(contentPane);
		
		okButton.addActionListener(this);
		cancelButton.addActionListener(this);
		
		this.pack();
		this.setLocationRelativeTo(null);  // centered on screen
		this.setVisible(true);
	}
	
	/**
	 * A method defined by the Actionlistener interface that allows this
	 * dialog to react accordingly when its buttons are clicked.
	 * 
	 * @param e An event heralding a click on one of this dialog's buttons.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();
		if (source == okButton || source == dataItemNameField) {
			String newName = dataItemNameField.getText();
			client.renameDataItem(data, newName);
			this.dispose();
		}
		if (source == cancelButton) {
			this.dispose();
		}
	}

	public void caretUpdate(CaretEvent e) {
		String dataSetName = dataItemNameField.getText();
		if (dataSetName.length() > 0) {
			okButton.setEnabled(true);
		} else {
			okButton.setEnabled(false);
		}
	}
}