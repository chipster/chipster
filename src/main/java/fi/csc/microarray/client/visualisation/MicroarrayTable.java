package fi.csc.microarray.client.visualisation;

import java.awt.Component;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.util.StringTokenizer;

import javax.swing.JComponent;
import javax.swing.JOptionPane;
import javax.swing.KeyStroke;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;
import javax.swing.table.TableCellRenderer;

import org.apache.log4j.Logger;
import org.jdesktop.swingx.JXTable;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.selection.RowChoiceEvent;
import fi.csc.microarray.databeans.DataBean;

/**
 * Table class which allows highlighting only the selected cell
 * not the entire row. Provides also copy-paste operations which
 * are compatible with most common spreadsheet editor like MS Excel 
 * and OpenOffice.org
 * 
 * @author mkoski
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
		KeyStroke copy = KeyStroke.getKeyStroke(KeyEvent.VK_C,ActionEvent.CTRL_MASK,false);
		KeyStroke paste = KeyStroke.getKeyStroke(KeyEvent.VK_V,ActionEvent.CTRL_MASK,false);
		this.registerKeyboardAction(this, "Copy", copy, JComponent.WHEN_FOCUSED);
		this.registerKeyboardAction(this, "Paste", paste, JComponent.WHEN_FOCUSED);
		
		systemClipboard = Toolkit.getDefaultToolkit().getSystemClipboard();				
		
		this.getSelectionModel().addListSelectionListener(new ListSelectionListener(){
			public void valueChanged(ListSelectionEvent e) {
				if(!doNotDispatchEvents){
					logger.debug("Row selection of the table changed");
					
					int[] selected = MicroarrayTable.this.getSelectedRows();
					int[] converted = new int[selected.length];
					
					for (int i = 0; i < selected.length; i++){
						converted[i] = MicroarrayTable.this.convertRowIndexToModel(selected[i]);
					}
					
					boolean tmp = doNotDispatchEvents;
					doNotDispatchEvents = true;
					application.getSelectionManager().getRowSelectionManager(MicroarrayTable.this.data).setSelected(
						converted, MicroarrayTable.this);
					
					doNotDispatchEvents = tmp;
				}
			}			
		});
		
		application.addPropertyChangeListener(this);				
	}
	
	/**
	 * Overrides prepareRenderer function so that the new CellRenderer
	 * changes only the background of the selected cell not the whole 
	 * row as the default CellRenderer does.
	 */
	@Override
	public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
	    
		// Returns cell under the event location
		Component cell = super.prepareRenderer(renderer, row, column);
	    
		// TODO The background coloring fails if selection mode is
		// set to MULTIPLY_INTERVAL_SELECTION mode. That's why the MULTIPLY_INTERVAL_SELECTION
		// mode is disabled. If the MULTIPLY_INTERVAL_SELECTION is needed in 
		// future this bug should be fixed.
		
	    if (this.isColumnSelected(column) && this.isRowSelected(row)){
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
	 * Prevents the use of MULTIPLE_INTERVAL_SELECTION. Before using multiple selection
	 * the cell background coloring problem occured in prepareRenderer method should be
	 * fixed
	 *
	@Override
	public void setSelectionMode(int selectionMode) {
		if(selectionMode == ListSelectionModel.MULTIPLE_INTERVAL_SELECTION){
			throw new IllegalArgumentException("MULTIPLE_INTERVAL_SELECTION cannot be used in MicroarrayTable");
		}
	}*/
	
	public void sendEvents(boolean value){
		doNotDispatchEvents = !value;
	}
	
	/**
	 * Copies selected cells to clipboard. Columns are separeted using 
	 * tab character and rows are separeted by using new line character. 
	 * This is the same format as the most common spreadsheet editor like 
	 * MS Excel or OpenOffice.org uses so copying is MS Excel compatible.
	 *
	 */
	public void copy(){
		StringBuffer selectedValues = new StringBuffer();

		int numColsSelected = this.getSelectedColumnCount();
		int numRowsSelected = this.getSelectedRowCount();
		
		int[] rowsSelected = this.getSelectedRows();
		int[] colsSelected = this.getSelectedColumns();

		for (int row = 0; row < numRowsSelected; row++)
		{
			for (int col = 0; col < numColsSelected; col++)
			{
				selectedValues.append(this.getValueAt(rowsSelected[row], colsSelected[col]));
				if (col < numColsSelected-1){
					selectedValues.append("\t");	// Adds tab to separate columns 
				}
			}
			selectedValues.append("\n");			// Adds line break to separate rows
		}
		StringSelection stringToClipboard = new StringSelection(selectedValues.toString());
		systemClipboard.setContents(stringToClipboard, stringToClipboard);
	}
	
	/**
	 * Pastes values to spreadsheet. Values must be in MS Excel format. 
	 * Columns are separeted by tab character and rows are separeted by 
	 * new line character. Does not edit the cell value if the cell is
	 * not editable.
	 */
	public void paste(){
		
		int startRow = (this.getSelectedRows())[0];
		int startCol = (this.getSelectedColumns())[0];
		
		try{
			String fromClipboard = (String)(systemClipboard.getData(DataFlavor.stringFlavor));
			
			StringTokenizer tokenizedValues = new StringTokenizer(fromClipboard,"\n");
			for(int row = 0; tokenizedValues.hasMoreTokens(); row++){
				String rowString = tokenizedValues.nextToken();
				StringTokenizer tokenizedRowValues = new StringTokenizer(rowString,"\t");
				
				for(int col = 0; tokenizedRowValues.hasMoreTokens(); col++){
					String cellValue = (String)tokenizedRowValues.nextToken();
					if ((startRow + row < this.getRowCount()) && (startCol + col < this.getColumnCount()) && this.isCellEditable(startRow + row, startCol + col)){
						this.setValueAt(cellValue, startRow + row, startCol + col);
					}
				}
			}
		} catch(IOException ioe){
			JOptionPane.showMessageDialog(this, "Error occured while retrieving data from clipboard", "Clipboard error", JOptionPane.ERROR_MESSAGE);
		} catch(UnsupportedFlavorException ufe){
			JOptionPane.showMessageDialog(this, "Requested data is not available", "Clipboard error", JOptionPane.ERROR_MESSAGE);
		} catch(IllegalStateException ise){
			JOptionPane.showMessageDialog(this, "Clipboard is currently unavailable", "Clipboard error", JOptionPane.ERROR_MESSAGE);
		}
	}

	public void actionPerformed(ActionEvent e) {
		if (e.getActionCommand().compareTo("Copy") == 0){
			copy();
		}
		if (e.getActionCommand().compareTo("Paste") == 0){
			paste();
		}
	}

	public void propertyChange(PropertyChangeEvent evt) {
		if(logger.isDebugEnabled()){
			if (evt instanceof RowChoiceEvent){
				logger.debug("Got a RowChoiceEvent from \n" + evt.getSource());
			}
		}
		
		if (evt instanceof RowChoiceEvent &&
				!(evt.getSource() == this) &&
				((RowChoiceEvent)evt).getData() == data){
			
			logger.debug("RowChoiceEvent not from Spreadsheet");
			
			updateSelectionsFromApplication();
		}
	}
	
	public void updateSelectionsFromApplication(){
		
		logger.debug("Updating selections from application");
		boolean tmp = doNotDispatchEvents;
		doNotDispatchEvents = true;
		this.getSelectionModel().clearSelection();
		/*int[][] intervals = application.getSelectionManager().getRowSelectionManager().getSelectedIntervals();		
		for (int i = 0; i < intervals.length; i++){
			this.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
			this.getSelectionModel().addSelectionInterval(intervals[i][0], intervals[i][1]);
		}*/
		this.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
		for (int row : application.getSelectionManager().getRowSelectionManager(data).getSelectedRows()){
			this.changeSelection(this.convertRowIndexToView(row), 0, true, false);
		}
		doNotDispatchEvents = tmp;
	}
}
