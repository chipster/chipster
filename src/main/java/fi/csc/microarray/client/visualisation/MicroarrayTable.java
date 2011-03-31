package fi.csc.microarray.client.visualisation;

import java.awt.Component;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableCellRenderer;

import org.apache.log4j.Logger;
import org.jdesktop.swingx.JXTable;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.selection.SelectionEvent;
import fi.csc.microarray.databeans.DataBean;

/**
 * Table class which allows highlighting only the selected cell not the entire
 * row. Provides also copy-paste operations which are compatible with most
 * common spreadsheet editor like MS Excel and OpenOffice.org
 * 
 * @author Mikko Koski, Aleksi Kallio
 * 
 */
public class MicroarrayTable extends JXTable implements ActionListener, PropertyChangeListener {

	private Clipboard systemClipboard;

	private static ClientApplication application = Session.getSession().getApplication();

	private static Logger logger = Logger.getLogger(MicroarrayTable.class);

	private boolean doNotDispatchEvents = true;
	private DataBean data;

	public MicroarrayTable(DataBean data) {

		this.data = data;

		// Register ctrl+c and ctrl+v
		KeyStroke copy = KeyStroke.getKeyStroke(KeyEvent.VK_C, ActionEvent.CTRL_MASK, false);
		KeyStroke paste = KeyStroke.getKeyStroke(KeyEvent.VK_V, ActionEvent.CTRL_MASK, false);
		this.registerKeyboardAction(this, "Copy", copy, JComponent.WHEN_FOCUSED);
		this.registerKeyboardAction(this, "Paste", paste, JComponent.WHEN_FOCUSED);

		systemClipboard = Toolkit.getDefaultToolkit().getSystemClipboard();

		this.getSelectionModel().addListSelectionListener(new ListSelectionListener() {
			public void valueChanged(ListSelectionEvent e) {
				if (!doNotDispatchEvents) {
					logger.debug("Row selection of the table changed");

					int[] selected = MicroarrayTable.this.getSelectedRows();
					int[] converted = new int[selected.length];

					for (int i = 0; i < selected.length; i++) {
						converted[i] = MicroarrayTable.this.convertRowIndexToModel(selected[i]);
					}

					boolean tmp = doNotDispatchEvents;
					doNotDispatchEvents = true;
					application.getSelectionManager().getRowSelectionManager(
							MicroarrayTable.this.data).setSelection(converted, MicroarrayTable.this);

					doNotDispatchEvents = tmp;
				}
			}
		});

		application.addClientEventListener(this);
	}

	/**
	 * Overrides prepareRenderer function so that the new CellRenderer changes
	 * only the background of the selected cell not the whole row as the default
	 * CellRenderer does.
	 */
	@Override
	public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {

		// Returns cell under the event location
		Component cell = super.prepareRenderer(renderer, row, column);

		// TODO The background coloring fails if selection mode is
		// set to MULTIPLY_INTERVAL_SELECTION mode. That's why the
		// MULTIPLY_INTERVAL_SELECTION
		// mode is disabled. If the MULTIPLY_INTERVAL_SELECTION is needed in
		// future this bug should be fixed.

		if (this.isColumnSelected(column) && this.isRowSelected(row)) {
			// Cell is selected
			cell.setBackground(this.getSelectionBackground());
			cell.setForeground(this.getSelectionForeground());
		} else {
			cell.setBackground(this.getBackground());
			cell.setForeground(this.getForeground());
		}
		return cell;
	}

	/**
	 * Prevents the use of MULTIPLE_INTERVAL_SELECTION. Before using multiple
	 * selection the cell background coloring problem occured in prepareRenderer
	 * method should be fixed
	 * 
	 * @Override public void setSelectionMode(int selectionMode) {
	 *           if(selectionMode ==
	 *           ListSelectionModel.MULTIPLE_INTERVAL_SELECTION){ throw new
	 *           IllegalArgumentException(
	 *           "MULTIPLE_INTERVAL_SELECTION cannot be used in MicroarrayTable"
	 *           ); } }
	 */

	public void sendEvents(boolean value) {
		doNotDispatchEvents = !value;
	}

	/**
	 * Copies selected cells to clipboard. Columns are separeted using tab
	 * character and rows are separeted by using new line character. This is the
	 * same format as the most common spreadsheet editor like MS Excel or
	 * OpenOffice.org uses so copying is MS Excel compatible.
	 * 
	 * Value of a single cell is copied without tab or new line.
	 * 
	 */
	public void copy() {
		StringBuffer selectedValues = new StringBuffer();

		int numColsSelected = this.getSelectedColumnCount();
		int numRowsSelected = this.getSelectedRowCount();

		int[] rowsSelected = this.getSelectedRows();
		int[] colsSelected = this.getSelectedColumns();

		// deal with a single cell separately to avoid adding unnecessary tab
		// and line break
		if (numColsSelected == 1 && numRowsSelected == 1) {
			selectedValues.append(this.getValueAt(rowsSelected[0], colsSelected[0]));
		}

		// multiple cells selected
		else {
			for (int row = 0; row < numRowsSelected; row++) {
				for (int col = 0; col < numColsSelected; col++) {
					selectedValues.append(this.getValueAt(rowsSelected[row], colsSelected[col]));
					if (col < numColsSelected - 1) {
						selectedValues.append("\t"); // Adds tab to separate
														// columns
					}
				}
				selectedValues.append("\n"); // Adds line break to separate rows
			}
		}
		StringSelection stringToClipboard = new StringSelection(selectedValues.toString());
		systemClipboard.setContents(stringToClipboard, stringToClipboard);
	}

	/**
	 * Pastes values to spreadsheet. Values must be in MS Excel format. Columns
	 * are separeted by tab character and rows are separeted by new line
	 * character. Does not edit the cell value if the cell is not editable.
	 * 
	 * Paste is done using model coordinates to avoid messing things up when
	 * pasting to sorted column.
	 * 
	 */
	public void paste() {

		int[] selectedRowsView = this.getSelectedRows();
		int[] selectedColsView = this.getSelectedColumns();

		// sanity check that target selection exists
		if (selectedRowsView.length < 1) {
			logger.warn("paste with no selected rows");
			return;
		}
		if (selectedColsView.length < 1) {
			logger.warn("paste with no selected columns");
		}


		// get data from clipboard
		List<List<String>> clipboardRows = new ArrayList<List<String>>();
		try {

			String fromClipboard = (String) (systemClipboard.getData(DataFlavor.stringFlavor));

			// use negative limit to prevent split from removing ignoring
			// sequential delimeters meaning several empty cells
			for (String row : fromClipboard.split("\n", -1)) {
				clipboardRows.add(Arrays.asList(row.split("\t", -1)));
			}
			// If there are more than one lines, split with -1 limit generates
			// one extra row in the end of the row. It's removed here.
			if (clipboardRows.size() > 1) {
				clipboardRows.remove(clipboardRows.size() - 1);
			}

		} catch (Exception e) {
			application.showDialog("Paste failed.", "Could not get the contents of the clipboard.", e.toString(), Severity.INFO, true);
			return;
		}

		
		// figure out the paste target cells

		// for selecting modified cells afterwards
		ArrayList<int[]> modifiedCells = new ArrayList<int[]>();

		// use model coords (to avoid sorting messing things during paste)
		int[] targetRowsModel;
		int[] targetColsModel;

		
		// Only one cell selected. The selected cell is the upper left corner of the
		// target area rectangle. The size of the rectange is limited either by the
		// contents of the clipboard or the size of the table, which ever happens to 
		// be smaller.
		if (selectedRowsView.length == 1 && selectedColsView.length == 1) {

			// figure out the height (rows) of the target area
			int targetAreaHeight;

			// clipboard contents limits rows
			if (this.getRowCount() - selectedRowsView[0] >= clipboardRows.size()) {
				targetAreaHeight = clipboardRows.size();
			} 

			// size of the table limits rows
			else {
				targetAreaHeight = this.getRowCount() - selectedRowsView[0];
			}


			// figure out the width (cols) of the target area
			int targetAreaWidth;

			// get the max width of the clipboard content
			int clipboardContentWidthMax = 0;
			for (List<String> clipboardRow : clipboardRows) {
				if (clipboardRow.size() > clipboardContentWidthMax) {
					clipboardContentWidthMax = clipboardRow.size();
				}
			}

			// clipboard contents limit the cols
			if (this.getColumnCount() - selectedColsView[0] >= clipboardContentWidthMax) {
				targetAreaWidth = clipboardContentWidthMax;
			} 

			// table size limits cols
			else {
				targetAreaWidth = this.getColumnCount() - selectedColsView[0];
			}

			// get the model coords for target area
			targetRowsModel = new int[targetAreaHeight];
			for (int i = 0; i < targetAreaHeight; i++) {
				targetRowsModel[i] = this.convertRowIndexToModel(selectedRowsView[0] + i);
			}
			targetColsModel = new int[targetAreaWidth];
			for (int i = 0; i < targetAreaWidth; i++) {
				targetColsModel[i] = this.convertColumnIndexToModel(selectedColsView[0] + i);
			}
		}

		// Multiple cells selected. The target cells are filled with the clipboard content, 
		// If the clipboard content is smaller than the target area, it's duplicated as
		// many times as needed.
		else {

			// get model coords for the target area
			targetRowsModel = new int[selectedRowsView.length];
			for (int i = 0; i < selectedRowsView.length; i++) {
				targetRowsModel[i] = this.convertRowIndexToModel(selectedRowsView[i]);
			}
			targetColsModel = new int[selectedColsView.length];
			for (int i = 0; i < selectedColsView.length; i++) {
				targetColsModel[i] = this.convertColumnIndexToModel(selectedColsView[i]);
			}
		}


		// fill the target area with clipboard content
		for (int i = 0; i < targetRowsModel.length; i++) {
			for (int j = 0; j < targetColsModel.length; j++) {
				List<String> clipboardRow = clipboardRows.get(i % clipboardRows.size());

				int targetRowModel = targetRowsModel[i];
				int targetColModel = targetColsModel[j];
				int targetRowView = this.convertRowIndexToView(targetRowModel);
				int targetColView = this.convertColumnIndexToView(targetColModel);

				if (this.isCellEditable(targetRowView, targetColView)) {
					this.setValueAt(clipboardRow.get(j % clipboardRow.size()), targetRowView,targetColView);
					modifiedCells.add(new int[] { targetRowModel, targetColModel });
				}
			}
		}

		// select the modified cells (though the ListSelectionModel doesn't work too well
		// when selecting multiple columns)
		this.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		this.clearSelection();
		for (int[] modifiedCell : modifiedCells) {
			int selectionRowView = this.convertRowIndexToView(modifiedCell[0]);
			int selectionColView = this.convertColumnIndexToView(modifiedCell[1]);
			this.changeSelection(selectionRowView, selectionColView, true, false);
		}
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().compareTo("Copy") == 0) {
			copy();
		}
		if (e.getActionCommand().compareTo("Paste") == 0) {
			paste();
		}
	}

	public void propertyChange(PropertyChangeEvent evt) {
		if (logger.isDebugEnabled()) {
			if (evt instanceof SelectionEvent) {
				logger.debug("Got a SelectionEvent from \n" + evt.getSource());
			}
		}

		if (evt instanceof SelectionEvent && !(evt.getSource() == this)
				&& ((SelectionEvent) evt).getData() == data) {

			logger.debug("SelectionEvent not from Spreadsheet");

			updateSelectionsFromApplication();
		}
	}

	public void updateSelectionsFromApplication() {

		logger.debug("Updating selections from application");
		boolean tmp = doNotDispatchEvents;
		doNotDispatchEvents = true;
		this.getSelectionModel().clearSelection();
		/*
		 * int[][] intervals =
		 * application.getSelectionManager().getRowSelectionManager
		 * ().getSelectedIntervals(); for (int i = 0; i < intervals.length;
		 * i++){
		 * this.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION
		 * ); this.getSelectionModel().addSelectionInterval(intervals[i][0],
		 * intervals[i][1]); }
		 */
		this.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		for (int row : application.getSelectionManager().getRowSelectionManager(data)
				.getSelectionAsRows()) {
			this.changeSelection(this.convertRowIndexToView(row), 0, true, false);
		}
		doNotDispatchEvents = tmp;
	}
}
