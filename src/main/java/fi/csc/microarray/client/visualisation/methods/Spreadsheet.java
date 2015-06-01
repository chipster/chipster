package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.table.DefaultTableModel;

import org.jdesktop.swingx.hyperlink.LinkModel;
import org.jdesktop.swingx.hyperlink.LinkModelAction;
import org.jdesktop.swingx.renderer.DefaultTableRenderer;
import org.jdesktop.swingx.renderer.HyperlinkProvider;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.selection.IntegratedSelectionManager;
import fi.csc.microarray.client.visualisation.ExtendedJXTable;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.RestrictModifier;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.basic.BasicModule;

/**
 * A GUI component for showing tabular data. Supports sorting and
 * other extended functions with underlying JXTable.
 * 
 * @author Janne Käki, Mikko Koski, Aleksi Kallio
 *
 */


public class Spreadsheet extends Visualisation {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}	

	/**
	 * PopupMenu for spread sheet view. Allows copy operation and annotating 
	 * using Bioconductor
	 *
	 */
	public class SpreadsheetPopupMenu extends JPopupMenu {

		public SpreadsheetPopupMenu(final ExtendedJXTable table, Module primaryModule) {

			JMenuItem copyMenuItem = new JMenuItem();
			copyMenuItem.setText("Copy");
			copyMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					table.copy();
				}
			});

			this.add(copyMenuItem);
			this.add(new JSeparator());


			JMenuItem filterMenuItem = new JMenuItem();
			filterMenuItem.setText("Create dataset from selected");
			filterMenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					VisualisationUtilities.filterBySelection(getFrame().getDatas(), false);
				}
			});
			this.add(filterMenuItem);
			
			JMenuItem filter2MenuItem = new JMenuItem();
			filter2MenuItem.setText("Create dataset from unselected");
			filter2MenuItem.addActionListener(new java.awt.event.ActionListener() {
				public void actionPerformed(java.awt.event.ActionEvent e) {
					VisualisationUtilities.filterBySelection(getFrame().getDatas(), true);
				}
			});
			this.add(filter2MenuItem);
		}
	}

	private final int COLUMNS_REQUIRES_SCROLLING = 8;

	private ExtendedJXTable table;

	/**
	 * Creates a new TablePanel, which is dataset specific.
	 * 
	 * @param datas The dataset for which to create the table.
	 * @throws Exception 
	 * @throws MicroarrayException 
	 */
	@Override
	public JComponent getVisualisation(final DataBean data) throws Exception {
		JPanel panel = new JPanel(new BorderLayout());				
		Module primaryModule = Session.getSession().getPrimaryModule();

		// Figure out column names
		String[] columnTitles;
		Object[][] rowData;
		List<Boolean> linkableFlags;
		int rowCount;
		int columnCount;
		QueryResult columnsFeature = data.queryFeatures("restrict(/column/*)");
		try (Table columns = columnsFeature.asTable()) {

			if (columns == null) {
				columnTitles = new String[] { "Info" };
				rowData = new String[][] { new String[] { DataManager.DATA_NA_INFOTEXT }};
				linkableFlags = new LinkedList<Boolean>();
				linkableFlags.add(false);
				rowCount = 1;
				columnCount = 1;

			} else {

				columnTitles = new String[columns.getColumnCount()];
				int counter = 0;
				for (String column : columns.getColumnNames()) {
					columnTitles[counter] = column;
					counter++;
				}
				columnCount = columns.getColumnCount();

				// Check which columns need hyperlinking
				linkableFlags = primaryModule.flagLinkableColumns(columns, data);

				// Count data rows
				try (Table rowCounter = data.queryFeatures("/column/*").asTable()) {
					rowCount = 0;
					while (rowCounter.nextRow()) {
						rowCount++;
					}
				}

				// Create actual tabular data
				rowData = new Object[RestrictModifier.RESTRICT_TO_ROWS < rowCount ? RestrictModifier.RESTRICT_TO_ROWS : rowCount][columns.getColumnCount()];
				int row = 0;
				while (columns.nextRow()) {
					int column = 0;
					for (String columnName : columns.getColumnNames()) {

						Object value = columns.getValue(columnName);
						ExtendedCellValue cell;

						IntegratedEntity linkedEntity = null;
						if (linkableFlags.get(column)) {
							// This cell value is linkable
							linkedEntity = primaryModule.createLinkableEntity(columns, data);
						}

						if (value instanceof Float) {
							cell = new ExtendedCellValue(columns.getStringValue(columnName), (Float)value, linkedEntity);

						} else {
							cell = new ExtendedCellValue(columns.getStringValue(columnName), null, linkedEntity);
						}
						rowData[row][column] = cell;
						column++;
					}
					row++;
				}
			}
		}

		// Create the table component
		table = new ExtendedJXTable(data);
		DefaultTableModel tableModel = new DefaultTableModel(rowData, columnTitles) {			
			@Override
			public boolean isCellEditable(int row, int column){
				return false;
			}			
		};
		table.setModel(tableModel);

		// Initialise support for hyperlinks, if needed
		for (int i = 0; i < linkableFlags.size(); i++) {

			if (linkableFlags.get(i)) {
				LinkModelAction<ExtendedCellValue> linkAction = new LinkModelAction<ExtendedCellValue>() {
					public void actionPerformed(ActionEvent e) {
						IntegratedSelectionManager selectionManager = application.getSelectionManager().getSelectionManager(data);
						selectionManager.setPointSelection(this.target.getLinkedEntity(), this);
					}
				};

				table.getColumn(i).setCellRenderer(new DefaultTableRenderer(new HyperlinkProvider(linkAction)));
			}
		}

		// Set look and feel aspects of the table 
		table.setColumnControlVisible(true);
		JScrollPane tableScroller = new JScrollPane(table);
		table.setBackground(java.awt.Color.white);
		table.setHorizontalScrollEnabled(columnCount > COLUMNS_REQUIRES_SCROLLING);

		// Initialise support for popups
		table.addMouseListener(new MouseAdapter(){
			@Override
			public void mousePressed(MouseEvent e) {
				maybeShowPopup(e);        
			}
			@Override
			public void mouseReleased(MouseEvent e) {
				maybeShowPopup(e);        
			}
			private void maybeShowPopup(MouseEvent e) {
				if (e.isPopupTrigger()) {
					JPopupMenu popup = new SpreadsheetPopupMenu(table, Session.getSession().getPrimaryModule());
					popup.show(e.getComponent(), e.getX(), e.getY());
				}
			}
		});

		// Make visible and activate
		panel.add(new JLabel("Showing " + rowData.length + " rows of " + rowCount), BorderLayout.NORTH);
		panel.add(tableScroller, BorderLayout.CENTER);
		table.updateSelectionsFromApplication();
		table.sendEvents(true);

		return panel; 
	}


	public static class ExtendedCellValue extends LinkModel implements Comparable {

		private String value;
		private Float numericValue;
		private IntegratedEntity linkedEntity;

		public ExtendedCellValue(String value, Float numericValue, IntegratedEntity linkedEntity) {
			this.value = value;
			this.numericValue = numericValue;
			this.linkedEntity = linkedEntity;
		}

		public IntegratedEntity getLinkedEntity() {
			return linkedEntity;
		}

		@Override
		public boolean equals(Object otherObj) {
			if (!(otherObj instanceof ExtendedCellValue)) {
				throw new IllegalArgumentException("cannot compare to instance of " + otherObj.getClass().getSimpleName());
			}

			ExtendedCellValue other = (ExtendedCellValue)otherObj;
			if (this.numericValue != null && other.numericValue != null) {
				return numericValue.equals(other.numericValue);

			} else {
				return value.equals(other.value);
			}

		}

		@Override
		public int hashCode() {
			return numericValue != null ? numericValue.hashCode() : value.hashCode();
		}

		@Override
		public String toString() {
			return getText();
		}

		@Override
		public String getText() {
			return value; // always use the String value for showing the value
		}

		@Override
		public int compareTo(Object otherObj) {
			if (!(otherObj instanceof ExtendedCellValue)) {
				throw new IllegalArgumentException("cannot compare to instance of " + otherObj.getClass().getSimpleName());
			}

			ExtendedCellValue other = (ExtendedCellValue)otherObj;
			if (this.numericValue != null && other.numericValue != null) {
				return numericValue.compareTo(other.numericValue);

			} else {
				return value.compareTo(other.value);
			}
		}

	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.hasTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES) || bean.hasTypeTag(BasicModule.TypeTags.TABLE_WITHOUT_COLUMN_NAMES);
	}

	@Override
	public void removeVisualisation(){
		application.removeClientEventListener(table);
	}
}
