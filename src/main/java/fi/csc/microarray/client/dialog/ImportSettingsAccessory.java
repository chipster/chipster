package fi.csc.microarray.client.dialog;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Set;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;

import fi.csc.microarray.client.dataimport.ImportUtils;

/**
 * Accessory component for the JFileChooser to allow direct input of dataset
 * folder.
 * 
 * @author Petri KlemelÃ¤
 * 
 */
public class ImportSettingsAccessory extends JPanel implements ActionListener {

	private static final String SKIP_TEXT = "Import directly if possible";
	private JCheckBox skipCheckBox;
	private JComboBox folderNameCombo;
	private String importFolder;

	public ImportSettingsAccessory(JFileChooser fileChooser) {

		this.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weightx = 1;
		
		c.gridx = 0;
		c.gridy = 0;
		c.insets.set(0, 5, 0, 0);
		JLabel folderNameLabel = new JLabel("Insert in folder");
		this.add(folderNameLabel, c);

		c.gridy++;
		Set<String> folderNameList = ImportUtils.getFolderNames(true);
		folderNameCombo = new JComboBox(folderNameList.toArray());
		folderNameCombo.setEditable(true);
		this.add(folderNameCombo, c);

		c.gridy++;
		c.insets.set(10, 5, 0, 0);
		c.weighty = 1;
		skipCheckBox = new JCheckBox(SKIP_TEXT);
		skipCheckBox.setSelected(true);
		this.add(skipCheckBox, c);
		
		// listen to file chooser for events that are related to this accessory
		fileChooser.addActionListener(this);
	}

	public boolean skipActionChooser() {
		return skipCheckBox.isSelected();
	}

	public void setDefaults() {
		skipCheckBox.setSelected(true);
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().equals(JFileChooser.APPROVE_SELECTION)) {
			this.setImportFolder(folderNameCombo.getSelectedItem().toString());
		}
	}

	public void setImportFolder(String folder) {
		this.importFolder = folder;
	}

	public String getImportFolder() {
		return this.importFolder;
	}
}