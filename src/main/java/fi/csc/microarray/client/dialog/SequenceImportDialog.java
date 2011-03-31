package fi.csc.microarray.client.dialog;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.ByteArrayInputStream;
import java.util.LinkedList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.databeans.DataBean;

@SuppressWarnings("serial")
public class SequenceImportDialog extends JDialog implements CaretListener, ActionListener {

	private static final String OPERATION_ID_FOR_MULTIPLE_IDENTIFIERS = "retrieve-sequences-from-database.sadl";
	private static final String OPERATION_ID_FOR_SINGLE_IDENTIFIER = "importseq.sadl";

	private static final String STYLE_FOR_HINTS = "style=\"font-size: 90%\"";
	
	private final Dimension BUTTON_SIZE = new Dimension(70, 25);

	private static Logger logger = Logger.getLogger(SequenceImportDialog.class);
	private ClientApplication application;

	private JLabel textLabel;
	private JTextField beginField;
	private JTextField endField;
	private JTextArea textArea;
	private JButton okButton;
	private JButton cancelButton;
	private JComboBox dbNameCombo;

	private enum Databases {
		EMBL("EMBL", "embl"), EMBL_NEW("EMBL New", "emblnew"), PDB("PDB", "pdb_seq"), UNIPROT_SWISS(
				"UniProt / SwissProt", "swiss"), UNIPROT_TREMBL("UniProt / TrEMBL", "trembl"), ENSEMBL("Ensembl",
				"ensembl");

		private String name;
		private String value;

		private Databases(String name, String value) {
			this.name = name;
			this.value = value;
		}

		public String getValue() {
			return value;
		}

		public String toString() {
			return name;
		}
	}

	public SequenceImportDialog(ClientApplication clientApplication) {
		super(Session.getSession().getFrames().getMainFrame(), true);

		this.application = clientApplication;
		this.setTitle("Import sequence");
		this.setModal(true);

		// Layout
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.WEST;
		c.gridx = 0;
		c.gridy = 0;
		c.gridwidth = 2;
		this.setLayout(new GridBagLayout());

		// Database
		dbNameCombo = new JComboBox(Databases.values());
		dbNameCombo.setPreferredSize(new Dimension(150, 20));
		dbNameCombo.setBackground(Color.WHITE);
		dbNameCombo.setSelectedItem(Databases.UNIPROT_SWISS);
		c.insets.set(10, 10, 5, 10);
		c.gridy++;
		this.add(new JLabel("Database"), c);
		c.insets.set(0, 10, 10, 10);
		c.gridy++;
		this.add(dbNameCombo, c);

		// Text label
		textLabel = new JLabel("Sequence identifier(s)");
		c.anchor = GridBagConstraints.WEST;
		c.insets.set(10, 10, 5, 10);
		c.gridy++;
		this.add(textLabel, c);

		// Text area
		textArea = new JTextArea();
		textArea.setFont(new Font(Font.MONOSPACED, Font.PLAIN, 12));
		textArea.setBorder(BorderFactory.createEmptyBorder());
		textArea.addCaretListener(this);
		textArea.setText("CASA1_HUMAN\nCASA1_MOUSE");
		textArea.selectAll();
		JScrollPane areaScrollPane = new JScrollPane(textArea);
		areaScrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		areaScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		areaScrollPane.setBorder(BorderFactory.createLineBorder(Color.GRAY));
		areaScrollPane.setPreferredSize(new Dimension(400, 90));
		c.insets.set(0, 10, 0, 10);
		c.gridy++;
		this.add(areaScrollPane, c);

//		c.gridy++;
//		c.anchor = GridBagConstraints.EAST;
//		c.gridwidth = 1;
//		c.gridx++;
//		this.add(new JLabel("<html><p " + STYLE_FOR_HINTS + ">separated by any whitespace, ';' or ','</p></html>"), c);
//		c.anchor = GridBagConstraints.WEST;
//		c.gridx = 0;
//		c.gridwidth = 2;
//		c.insets.set(0, 10, 10, 10);


		// range panel
		beginField = new JTextField(4);
		endField = new JTextField(4);
		JPanel rangePanel = new JPanel(new GridBagLayout());
		GridBagConstraints rangeConstraints = new GridBagConstraints();
		rangeConstraints.weightx = 1.0;
		rangeConstraints.weighty = 1.0;
		rangeConstraints.anchor = GridBagConstraints.WEST;
		rangeConstraints.insets.set(0, 0, 0, 10);
		rangeConstraints.gridx = 0;
		rangeConstraints.gridy = 0;
		rangePanel.add(new JLabel("Start"), rangeConstraints);
		rangeConstraints.gridx++;
		rangePanel.add(beginField, rangeConstraints);
		rangeConstraints.gridx++;
		rangePanel.add(new JLabel("End"), rangeConstraints);
		rangeConstraints.gridx++;
		rangePanel.add(endField, rangeConstraints);
		rangeConstraints.gridx = 0;
		rangeConstraints.gridy++;
		rangeConstraints.gridwidth = 4;
		rangeConstraints.anchor = GridBagConstraints.EAST;
//		rangePanel.add(new JLabel("<html><p " + STYLE_FOR_HINTS + ">only applies for a single identifier</p></html>"), rangeConstraints);
		
		c.insets.set(10, 10, 0, 10);
		c.gridy++;
		c.fill = GridBagConstraints.NONE;
		this.add(rangePanel, c);
		
		c.insets.set(0, 10, 0, 10);
		c.gridy++;
		c.fill = GridBagConstraints.NONE;
		this.add(new JLabel("<html><p " + STYLE_FOR_HINTS + ">Only applies for a single identifier</p></html>"), c);
		
		// OK button
		okButton = new JButton("OK");
		okButton.setPreferredSize(BUTTON_SIZE);
		okButton.addActionListener(this);

		// Cancel button
		cancelButton = new JButton("Cancel");
		cancelButton.setPreferredSize(BUTTON_SIZE);
		cancelButton.addActionListener(this);

		// Buttons panel
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

		c.insets.set(10, 10, 10, 10);
		c.anchor = GridBagConstraints.SOUTHEAST;
		c.gridy++;
		c.fill = GridBagConstraints.NONE;

		this.add(keepButtonsRightPanel, c);

		// Show
		this.pack();
		this.setResizable(false);
		Session.getSession().getFrames().setLocationRelativeToMainFrame(this);

		// Default focus
		textArea.requestFocusInWindow();
		this.setVisible(true);
	}

	public void caretUpdate(CaretEvent e) {

	}

	public void actionPerformed(ActionEvent e) {

		if (e.getSource() == okButton) {

			// db
			String db = ((Databases) this.dbNameCombo.getSelectedItem()).getValue();

			// identifiers are separated by any number of whitespace ; or ,
			// also make sure that no empty strings are added as identifiers
			String identifiersText = this.textArea.getText().trim();
			List<String> idList = new LinkedList<String>();
			for (String id: identifiersText.split("[\\s|;|,]+")) {
				if (!id.isEmpty()) {
					idList.add(id);
				}
			}

			// range
			Integer sbegin;
			try {
				sbegin = Integer.parseInt(beginField.getText());
			} catch (NumberFormatException nfe) {
				sbegin = null;
			}
			Integer send;
			try {
				send = Integer.parseInt(endField.getText());
			} catch (NumberFormatException nfe) {
				send = null;
			}

			// run
			if (idList.size() > 0) {
				runImport(db, idList.toArray(new String[] {}), sbegin, send);
				this.dispose();
			}
		}

		// cancel
		else if (e.getSource() == cancelButton) {
			this.dispose();
		}
	}

	/**
	 * Run import job.
	 * 
	 * @param fileName
	 * @param folderName
	 * @param db
	 * @param id
	 */
	private void runImport(String db, String[] identifiers, Integer sbegin, Integer send) {
		try {
			// Create importing job
			logger.info("Importing sequence...");

			String identifiersString = "";
			for (String identifier : identifiers) {
				identifiersString += identifier + "\n";
			}

			// create the right operation, depending on if there are one ore more identifiers
			Operation operation;
			if (identifiers.length == 1) {
				operation = new Operation(application.getOperationDefinition(OPERATION_ID_FOR_SINGLE_IDENTIFIER), new DataBean[] {});
				operation.setParameter("sequence", db + ":" + identifiers[0]);
				operation.setParameter("sbegin", sbegin);
				operation.setParameter("send", send);


			} else {
				operation = new Operation(application.getOperationDefinition(OPERATION_ID_FOR_MULTIPLE_IDENTIFIERS), new DataBean[] {});
				operation.setParameter("source", db);
				operation.bindInputs(new DataBean[] { Session.getSession().getDataManager().createDataBean(
						"identifiers.txt", new ByteArrayInputStream(identifiersString.getBytes())) });
			}
			
			// run the task
			application.executeOperation(operation);

		} catch (Exception exc) {
			application.reportException(exc);
		}
	}
}
