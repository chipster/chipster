package fi.csc.microarray.client.dataimport;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;

import org.jdesktop.swingx.JXTable;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.dataimport.ImportItem.Action;
import fi.csc.microarray.client.dialog.ChipsterDialog;
import fi.csc.microarray.client.dialog.DialogInfo;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.DialogCloseListener;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.dialog.DialogInfo.Type;
import fi.csc.microarray.constants.VisualConstants;

public class ActionChooserScreen implements ActionListener, DialogCloseListener {

	private class GreyTableCellRenderer extends DefaultTableCellRenderer {
		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			setEnabled(false);
			return super.getTableCellRendererComponent(table, value, isSelected, hasFocus, row, column);
		}
	}

	public class ComboBoxRenderer extends JComboBox implements TableCellRenderer {

		public ComboBoxRenderer(ImportItem.Action[] items) {
			super(items);

		}

		public Component getTableCellRendererComponent(JTable table, Object value, boolean isSelected, boolean hasFocus, int row, int column) {
			if (isSelected) {
				setForeground(table.getSelectionForeground());
				super.setBackground(table.getSelectionBackground());
			} else {
				setForeground(table.getForeground());
				setBackground(table.getBackground());
			}

			// Select the current value
			setSelectedItem(value);
			return this;
		}
	}

	private static final String TITLE_TEXT = "Choose how to proceed with each file";
	private static final String INFO_TEXT = "A file can be imported directly as it is or you can use the Import tool to define the contents of the file. You can also decide not to import a file at all.";
	
	private static final int NAME_COLUMN_INDEX = 0;
	private static final int TYPE_COLUMN_INDEX = 1;
	private static final int ACTION_COLUMN_INDEX = 2;

	private JDialog dialog = new JDialog(Session.getSession().getFrames().getMainFrame(), "Import");

	private JCheckBox sameSettingsCheckBox;

	private JButton okButton;
	private JButton cancelButton;

	private JXTable table;
	private DefaultTableModel tableModel;

	private ImportSession importSession;

	public ActionChooserScreen(ImportSession importSession) {
		SwingClientApplication.setPlastic3DLookAndFeel(dialog);
		dialog.setPreferredSize(new Dimension(640, 480));
		dialog.setLocationByPlatform(true);

		this.importSession = importSession;

		table = this.getTable();
		JScrollPane scroll = new JScrollPane(table);

		// upper panel
		JLabel titleLabel = new JLabel("<html><p style=" + VisualConstants.HTML_DIALOG_TITLE_STYLE + ">" + TITLE_TEXT + "</p></html>", JLabel.LEFT);
		JLabel descriptionLabel = new JLabel("<html><p>" + INFO_TEXT + "</p></html>", JLabel.LEFT);
		GridBagConstraints c = new GridBagConstraints();
		JPanel upperPanel = new JPanel(new GridBagLayout());
		c.weightx = 1.0;
		c.weighty = 1.0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.insets.set(10, 10, 5, 10);
		c.gridx = 0;
		c.gridy = 0;
		upperPanel.add(titleLabel, c);
		c.gridy++;
		upperPanel.add(descriptionLabel, c);
		

		// lower panel
		JPanel buttonPanel = new JPanel();
		okButton = new JButton("  OK  ");
		cancelButton = new JButton("Cancel");

		okButton.addActionListener(this);
		cancelButton.addActionListener(this);

		JLabel sameSettingsLabel = new JLabel("<html><p>" + "When using Import tool to import more than one files, only define the contents of the first file and then apply the same settings for the rest of the files." + "</p></html>", JLabel.LEFT);
		sameSettingsLabel.setVerticalTextPosition(JLabel.TOP);
		//sameSettingsLabel.setPreferredSize(new Dimension(550, 40));

		sameSettingsCheckBox = new JCheckBox("Define file structure once and apply the same settings to all files");
		sameSettingsCheckBox.setEnabled(true);
		sameSettingsCheckBox.setSelected(true);
		sameSettingsCheckBox.setPreferredSize(new Dimension(550, 40));

		buttonPanel.setLayout(new GridBagLayout());
		GridBagConstraints g = new GridBagConstraints();
		g.anchor = GridBagConstraints.NORTHWEST;
		g.gridx = 0;
		g.gridy = 0;
		g.weightx = 0.0;

		g.insets = new Insets(5, 5, 10, 5);
		buttonPanel.add(sameSettingsCheckBox, g);
		g.insets = new Insets(5, 0, 10, 5);
		g.gridy++;
		g.anchor = GridBagConstraints.EAST;
		buttonPanel.add(cancelButton, g);
		g.gridx++;
		buttonPanel.add(okButton, g);

		dialog.setLayout(new BorderLayout());
		dialog.add(upperPanel, BorderLayout.NORTH);
		dialog.add(scroll, BorderLayout.CENTER);
		dialog.add(buttonPanel, BorderLayout.SOUTH);
		dialog.pack();
		dialog.pack();
		dialog.setVisible(true);
	}

	/**
	 * Gets the table if it exists or creates a new one.
	 * 
	 * @return
	 */
	private JXTable getTable() {
		if (table == null) {

			class ActionChooserTableModel extends DefaultTableModel {

				public final String[] COLUMNS = new String[] { "Filename", "Detected type", "Action" };

				public ActionChooserTableModel(ImportSession importSession) {

					// populate the table
					Object[][] data = new Object[importSession.getItemCount()][];

					for (int row = 0; row < data.length; row++) {
						ImportItem item = importSession.getItemAtIndex(row);
						Object[] rowData = { item.getOutput().getName(), item.getType(), item.getAction() };
						data[row] = rowData;
					}

					setDataVector(data, COLUMNS);
				}

				@Override
				public boolean isCellEditable(int row, int column) {
					return column == NAME_COLUMN_INDEX || column == ACTION_COLUMN_INDEX;
				}

				@Override
				public int getRowCount() {
					return importSession.getItemCount();
				}

				@Override
				public Object getValueAt(int row, int column) {

					// Too big or negative row number
					if (row >= getRowCount() || row < 0) {
						throw new ArrayIndexOutOfBoundsException(row);
					}

					ImportItem item = importSession.getItemAtIndex(row);

					// Filename
					if (column == 0) {
						return item.getOutput().getName();
					}

					// Type
					else if (column == 1) {
						return item.getType().getDescription();
					}

					// Action
					else if (column == 2) {
						return item.getAction();
					}

					// Too big or negative column number
					else {
						throw new ArrayIndexOutOfBoundsException(row);
					}
				}

				@Override
				public void setValueAt(Object value, int row, int column) {
					ImportItem item = importSession.getItemAtIndex(row);

					// change action
					if (column == ACTION_COLUMN_INDEX) {
						item.setAction((Action) value);
					}

					// change filename
					else if (column == NAME_COLUMN_INDEX) {
						item.setFilename((String) value);
					}

					// others are illegal
					else {
						throw new IllegalArgumentException("Illegal column: " + column);
					}

					// is this needed?
					fireTableCellUpdated(row, column);

					// if (column != 0) {
					// super.setValueAt(aValue, row, column);
					// } else {
					// ImportItem item = importSession.getItemAtIndex(row);
					// item.setAction((Action)aValue);
					// }
				}
			}

			tableModel = new ActionChooserTableModel(importSession);

			table = new JXTable(tableModel);

			table.getColumnModel().getColumn(TYPE_COLUMN_INDEX).setCellRenderer(new GreyTableCellRenderer());

			table.getColumnModel().getColumn(ACTION_COLUMN_INDEX).setCellRenderer(new ComboBoxRenderer(ImportItem.Action.values()));
			table.getColumnModel().getColumn(ACTION_COLUMN_INDEX).setCellEditor(new DefaultCellEditor(new JComboBox(ImportItem.Action.values())));

			return table;
		} else {
			return table;
		}
	}


	public void actionPerformed(ActionEvent e) {

		if (e.getSource() == okButton) {

			// check for bad settings
			String warnings = "";
			for (int row = 0; row < importSession.getItemCount(); row++) {
				int index = table.convertRowIndexToModel(row);
				ImportItem item = importSession.getItemAtIndex(index);
				File file = item.getInput();
				if (file.getName().endsWith(".zip") && item.getAction().equals(ImportItem.Action.CUSTOM)) {
					warnings = warnings + "File " + file.getName() + " is a ZIP file, and cannot be handled by import tool.\n";
				}
			}

			if (warnings.length() > 0) {
				DialogInfo info = new DialogInfo(Severity.WARNING, "Attempting import of binary files", "Import tool can describe only files that contain text formatted tables. For example ZIP archived files cannot be used. See details for more information. Select Ok to proceed or Cancel to go back.", warnings, Type.OPTION); 
				ChipsterDialog.showDialog(null, info, DetailsVisibility.DETAILS_VISIBLE, false, this, null);
				return; // return to wait for users decision
			}

			doImports();
			
		} else if (e.getSource() == cancelButton) {
			// Close the frame
			dialog.dispose();
		}
	}

	public void dialogClosed(boolean okSelected) {
		if (okSelected) {
			doImports();
		}
	}

	private void doImports() {
		
		List<ImportItem> directImportDatas = new ArrayList<ImportItem>();
		
		// do imports
		for (int row = 0; row < importSession.getItemCount(); row++) {
			// These go very fast, but maybe it helps if there is a problem
			table.setRowSelectionInterval(row, row);

			int index = table.convertRowIndexToModel(row);
			ImportItem item = importSession.getItemAtIndex(index);

			if (item.getAction().equals(ImportItem.Action.DIRECT)) {

				directImportDatas.add(item);				

			} else {
				// Files with action IGNORE will be skipped and CUSTOM
				// files are imported a bit later
			}
		}
		
		Session.getSession().getApplication().importGroup(directImportDatas, importSession.getDestinationFolder());

		// Check if there were any custom files
		if (importSession.hasCustomFiles()) {
			importSession.setUseSameDescriptions(sameSettingsCheckBox.isSelected());
			((SwingClientApplication) Session.getSession().getApplication()).openImportTool(importSession);
		}

		// Close the frame
		dialog.dispose();
	}
}
