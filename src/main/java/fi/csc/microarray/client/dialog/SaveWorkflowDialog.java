package fi.csc.microarray.client.dialog;

import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import fi.csc.microarray.databeans.DataBean;

// FIXME remove this, it is not used
public class SaveWorkflowDialog extends JDialog
							    implements ActionListener, CaretListener {

	private static int instanceCounter = 0;
	
	private final Dimension BUTTON_SIZE = new Dimension(75, 25);
	
	private JTextField workflowNameField;
	private JButton okButton;
	private JButton cancelButton;
	
	/**
	 * Creates a new Save Workflow creation dialog and shows it.
	 * 
	 * @param owner The owning frame component of the dialog to be created.
	 * 				In this case, it would be the GUI mainframe class.
	 * @param data The dataset from thich the new workflow is derived.
	 */
	public SaveWorkflowDialog(DataBean data, Component initialLocation) {
		super((JFrame)null, "Save Workflow", true);
		
		workflowNameField = new JTextField(15);
		workflowNameField.setText("Workflow " + instanceCounter+1);
		
		okButton = new JButton("OK");
		okButton.setPreferredSize(BUTTON_SIZE);
		cancelButton = new JButton("Cancel");
		cancelButton.setPreferredSize(BUTTON_SIZE);
		
		JPanel contentPane = new JPanel(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.WEST;
		c.insets.set(4, 12, 4, 10);
		c.gridx = 0; c.gridy = 0;
		c.insets.set(4, 12, 4, 10);
		JLabel dataSetNameLabel = new JLabel("Workflow Name");
		contentPane.add(dataSetNameLabel, c);
		c.gridx = 1;
		c.gridwidth = 3;
		contentPane.add(workflowNameField, c);
		c.gridy = 1;
		c.gridwidth = 1;
		contentPane.add(okButton, c);
		c.gridx = 3;
		contentPane.add(cancelButton, c);
		this.setContentPane(contentPane);
		
		okButton.addActionListener(this);
		cancelButton.addActionListener(this);
		
		this.pack();
        this.setLocationRelativeTo(initialLocation);
		this.setVisible(true);
	}
	
	/**
	 * A method defined by the ActionListener interface that allows this
	 * dialog to react accordingly when its buttons are clicked.
	 * 
	 * @param e An event heralding a click on one of this dialog's buttons.
	 */
	public void actionPerformed(ActionEvent e) {
		throw new UnsupportedOperationException();
	}

	/**
	 * A method defined by the CaretListener interface that allows this
	 * dialog to react accordingly when the user moves the caret
	 * (writing point) in a text component.
	 * 
	 * @param e An event heralding a caret positiion change.
	 */
	public void caretUpdate(CaretEvent e) {
		String workflowName = workflowNameField.getText();
		if (workflowName.length() > 0) {
			okButton.setEnabled(true);
		} else {
			okButton.setEnabled(false);
		}
	}
}