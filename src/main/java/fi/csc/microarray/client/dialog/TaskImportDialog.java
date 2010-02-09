package fi.csc.microarray.client.dialog;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.parameter.ImportParameterPanel;
import fi.csc.microarray.exception.MicroarrayException;


public class TaskImportDialog extends JDialog implements ActionListener {

	private final Dimension BUTTON_SIZE = new Dimension(70, 25);

	private JLabel titleLabel;
	private JLabel descriptionLabel;
	private JLabel noteLabel;
	private JButton okButton;
	private JButton cancelButton;
	private JComboBox folderNameCombo;
	private ClientApplication application;
	private Operation operation;

	public TaskImportDialog(SwingClientApplication application, String databaseName, Operation operation) throws MicroarrayException {
		super(application.getMainFrame(), true);

		this.application = application;
		this.operation = operation;
		this.setTitle("Import");
		this.setModal(true);
		this.setPreferredSize(new Dimension(500, 300));

		// initialise components
		titleLabel = new JLabel("<html><p style=\"font-weight:bold;font-size:120%\">Import data from " + databaseName + "</p></html>");
		descriptionLabel = new JLabel("<html>" + operation.getDescription() + "</html>");
		noteLabel = new JLabel("<html><p style=\"font-style:italic\">It may take a while for the import task to finish.");

		folderNameCombo = new JComboBox(ImportUtils.getFolderNames(false).toArray());
		folderNameCombo.setEditable(true);

		okButton = new JButton("Import");
		okButton.setPreferredSize(BUTTON_SIZE);
		okButton.addActionListener(this);

		cancelButton = new JButton("Cancel");
		cancelButton.setPreferredSize(BUTTON_SIZE);
		cancelButton.addActionListener(this);

		ImportParameterPanel parameterPanel = new ImportParameterPanel(operation, null);

		JPanel keepButtonsRightPanel = new JPanel(new GridBagLayout());
		GridBagConstraints buttonConstraints = new GridBagConstraints();
		buttonConstraints.weightx = 1.0;
		buttonConstraints.weighty = 1.0;
		buttonConstraints.anchor = GridBagConstraints.EAST;
		buttonConstraints.insets.set(0, 0, 0, 8);
		keepButtonsRightPanel.add(cancelButton, buttonConstraints);
		buttonConstraints.gridx = GridBagConstraints.RELATIVE;
		buttonConstraints.insets.set(0, 0, 0, 0);
		keepButtonsRightPanel.add(okButton, buttonConstraints);
		
		
		// layout
		GridBagConstraints c = new GridBagConstraints();
		this.setLayout(new GridBagLayout());

		// title label
		c.weightx = 1.0;
		c.weighty = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.insets.set(10, 10, 5, 10);
		c.gridx = 0;
		c.gridy = 0;
		this.add(titleLabel, c);

		// description label
		c.gridy++;
		this.add(descriptionLabel, c);
		
		// parameter panel
		c.gridy++;
		c.weighty = 120.0;
		c.fill = GridBagConstraints.BOTH;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.insets.set(0, 10, 10, 10);
		this.add(parameterPanel, c);

		// note
		c.gridy++;
		c.weightx = 1.0;
		c.weighty = 1.0;
		c.insets.set(0, 10, 10, 10);
		c.fill = GridBagConstraints.HORIZONTAL;
		this.add(noteLabel, c);
		
		// buttons
		c.insets.set(10, 10, 10, 10);
		c.anchor = GridBagConstraints.SOUTHEAST;
		c.gridy++;
		c.fill = GridBagConstraints.NONE;
		this.add(keepButtonsRightPanel, c);


		// make visible
		this.pack();
		this.setLocationRelativeTo(application.getMainFrame());
		this.setVisible(true);
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
}
