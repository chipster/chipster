package fi.csc.microarray.client.dataimport;

import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.dataimport.events.ChipCountChangeEvent;
import fi.csc.microarray.client.dataimport.events.ChipNumberChangedEvent;
import fi.csc.microarray.client.dataimport.events.ColumnTypeChangeListener;
import fi.csc.microarray.client.dataimport.events.ColumnTypeChangeSupport;
import fi.csc.microarray.client.dataimport.events.ColumnTypeChangedEvent;

/**
 * Class which manages the column types and chip numbers.
 * 
 * The class also implements <code>ColumnTypeChangeSupport</code> class 
 * so it can be listened and GUI components can be updated when event is 
 * occured.
 * 
 * @author mkoski
 *
 */
public class ColumnTypeManager implements ColumnTypeChangeSupport {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(ColumnTypeManager.class);
	
	/**
	 * List of columns. The list is actually one bigger than the data column count
	 * because of the first column which represents line numbers. So after writing 
	 * data to list the first column should be ignored.
	 */
	private List<DataColumn> columns;
	
	/**
	 * Chips count
	 */
	private int chipCount;
	
	/**
	 * Column count
	 */
	private int columnCount;
	
	/**
	 * Columns before pattern selection. This allows undoing 
	 * if the pattern selection gives an unwanted result
	 */
	private List<DataColumn> columnsBeforePatternSelection;
	
	/**
	 * Listeners 
	 */
	private List<ColumnTypeChangeListener> columnTypeListeners;
	
	/**
	 * Gives a new column type manager with a given column count
	 * 
	 * @param columnCount column count
	 */
	public ColumnTypeManager(int columnCount){
		this.columnCount = columnCount;
		this.chipCount = 0;
		this.columns = new ArrayList<DataColumn>();
		this.columnTypeListeners = new ArrayList<ColumnTypeChangeListener>();
		this.columnsBeforePatternSelection = new ArrayList<DataColumn>();
	}
	
	/**
	 * Sets column count. Also updates the columnTypes and columnChipNumbers
	 * lists if needed.
	 * 
	 * @param columnCount
	 */
	public void setColumnCount(int columnCount){
		this.columnCount = columnCount;
		
		// Fill the list
		int columnIndex = 0;
		while(columns.size() <= columnCount){
			if(columnIndex == 0){
				// First column is the row number
				columns.add(DataColumn.createRowNumberColumn(columnIndex));
			}
			columns.add(DataColumn.createEmptyDataColumn(columnIndex));
			columnIndex++;
		}
	}

	/**
	 * Sets chip count
	 * 
	 * @param chipCount chip count
	 */
	public void setChipCount(int chipCount){
		this.chipCount = chipCount;
		fireChipCountChangeEvent(new ChipCountChangeEvent(this, chipCount));
	}
	
	/**
	 * Gets count of how many given type columns is selected from the table.
	 * 
	 * @param askedType
	 * @return count of how many given type columns is selected from the table
	 */
	public int getCountOfType(ColumnType askedType){
		int count = 0;
		for(DataColumn columnFromList : this.columns){
			if(columnFromList.getColumnType().equals(askedType)){
				count++;
			}
		}
		return count;
	}
	
	/**
	 * Returns the column names of sample columns. The size of the returned list
	 * is the same as return value of getCountOfType(ColumnType.SAMPLE_LABEL).
	 * 
	 * @return
	 */
	public List<String> getOriginalChipNames(){
		List<String> names = new ArrayList<String>();
		
		for(DataColumn column : this.columns){
			if(column.getColumnType().equals(ColumnType.SAMPLE_LABEL)){
				names.add(column.getOriginalName());
			}
		}
		return names;
	}
			
	
	/**
	 * Sets column type to given column index and removes the chip number.
	 * 
	 * @param  columnIndex column index number
	 * @param  type column type
	 * @param originalName 
	 * @throws IllegalArgumentException if user tries to change the type of the row number column 
	 * 		   or user tries to set a column type ROW_NUBMER
	 */
	public void setColumnType(int columnIndex, ColumnType type, String originalName){
		if(columnIndex == 0){
			throw new IllegalArgumentException("Can't change the type of the first column");
		} else if(ColumnType.ROW_NUMBER.equals(type)){
			throw new IllegalArgumentException("Can't change the column type to ROW_NUMBER");
		}
		
		// Sets column type and removes the chip number
		columns.set(columnIndex, new DataColumn(type, -1, columnIndex, originalName));
		
		// Every chip has one Sample. The count of Samples is equal to count of chips
		setChipCount(getCountOfType(ColumnType.SAMPLE_LABEL));
		
		// Notify listeners
		this.fireColumnTypeChangeEvent(new ColumnTypeChangedEvent(this, type, columnIndex));
	}
	
	/**
	 * Multiple column types can be selected using <code>ColumnTypePattern</code>.
	 * This is usually done when the "guess the rest" functionality is used. 
	 * The current column types are saved to a list and they can be returned 
	 * with the help of <code>undoPatternSelection</code>-method if the 
	 * result of the pattern selection is unwanted.
	 * 
	 * @param pattern pattern of the columns to be filled
	 * @param originalNames 
	 * @param ignore the unused columns at the beginning of the table
	 */
	public void selectColumnsFromPattern(ColumnTypePattern pattern, boolean ignoreFirstUnusedColumns, String[] originalNames){
		
		// Save the current column types
		this.columnsBeforePatternSelection.clear();
		this.columnsBeforePatternSelection.addAll(getColumns());
		
		boolean firstUsedColumnFound;
		
		if(ignoreFirstUnusedColumns){
			// Ignore first unused columns and start filling after those
			firstUsedColumnFound = false;
		} else {
			// Start filling from beginning of the table
			firstUsedColumnFound = true;
		}
		
		for(int columnNumber = 1; columnNumber < getColumnCount(); columnNumber++){
			if(!getColumnType(columnNumber).equals(ColumnType.UNUSED_LABEL)){
				firstUsedColumnFound = true;
			}
			if(firstUsedColumnFound){
				ColumnType type = pattern.getColumnTypeForIndex(columnNumber);
				if(type != null){
					
					String title = originalNames != null && originalNames.length >= columnNumber ? 
							originalNames[columnNumber] : null;
							
					setColumnType(columnNumber, type, title);
					setColumnChipNumber(columnNumber, getNextChipNumber(type));
				}
			}
		}
	}
	
	
	public void setColumns(List<DataColumn> newColumns){
		for(int columnNumber = 0; columnNumber < getColumnCount(); columnNumber++){
			if(columnNumber == 0){
				// Row number column. The row number is not changed, so no need to do anything
				continue;
			} else {
				DataColumn column = newColumns.get(columnNumber);
				setColumnType(columnNumber, column.getColumnType(), column.getOriginalName());
				setColumnChipNumber(columnNumber, column.getChipNumber());
			}
		}
	}
	
	/**
	 * Undo pattern selection.
	 */
	public void undoPatternSelection(){
		if(columnsBeforePatternSelection.isEmpty()){
			throw new IllegalArgumentException("Nothing to undo");
		}
		
		this.setColumns(columnsBeforePatternSelection);
		columnsBeforePatternSelection.clear();
	}
	
	/**
	 * Gets next chip number for type. If user has set numbers 1, 2 and 3, 
	 * the result is obviously 4 but if user has set chips 1, 3 and 4 the 
	 * result is 2
	 * 
	 * @param type
	 * @return next chip number for type.
	 */
	public int getNextChipNumber(ColumnType type){
		
		// Gets list of chip numbers which are in use
		List<Integer> chipNumbersInUse = new ArrayList<Integer>();
		for(int i = 0; i < columns.size(); i++){
			
			// Is the type same as asked type
			if(columns.get(i).getColumnType().equals(type)){
				
				// Yes it is! Add the chip number to the list
				chipNumbersInUse.add(columns.get(i).getChipNumber());
			}
		}
		
		logger.debug("Count of chip numbers set to " + type.toString() + " is " + chipNumbersInUse.size());
		
		// Start from 1 and search the first unused number
		for(int i = 1; i < columnCount; i++){
			if(!chipNumbersInUse.contains(i)){
				return i;
			}
		}
		
		// The first unused number found at the end of the list
		return columnCount;
	}
	
	/**
	 * Returns column type. If column type is not set for given column 
	 * return ColumnType.UNUSED_LABEL
	 * 
	 * @param columnIndex column number
	 * @return column type or UNUSED_LABEL if type is not set
	 */
	public ColumnType getColumnType(int columnIndex){
		return columns.get(columnIndex).getColumnType();
	}
	
	/**
	 * Gets columns chip number. If chip number is not set, negative value
	 * is returned
	 * 
	 * @param columnIndex
	 * @return columns chip number or negative value
	 */
	public int getColumnChipNumber(int columnIndex){
		return columns.get(columnIndex).getChipNumber();
	}
	
	/**
	 * Gets column count. Includes the first row number column.
	 * 
	 * @return column count
	 */
	public int getColumnCount(){
		return columnCount;
	}
	
	/**
	 * Chip count
	 * 
	 * @return chip count
	 */
	public int getChipCount(){
		return chipCount;
	}

	/**
	 * Short description of column type counts.
	 */
	public String toString(){
		StringBuffer string = new StringBuffer();
		string.append("ColumnTypeManager [column count: "+ getColumnCount() +" UNUSED: " + getCountOfType(ColumnType.UNUSED_LABEL) + " IDENTIFIERS: " + getCountOfType(ColumnType.IDENTIFIER_LABEL) + " SAMPLES: " + getCountOfType(ColumnType.SAMPLE_LABEL) + "]");
		return string.toString();
	}
	
	/**
	 * Has user selected less columns for certain type than chips exists?
	 * 
	 * @param type
	 * @return
	 */
	public boolean isLessTypesSelectedThanChips(ColumnType type){
		return getCountOfType(type) < getChipCount();
	}
	
	/**
	 * Has user selected as many columns for certain type than chips exists?
	 * 
	 * @param type
	 * @return
	 */
	public boolean isAsManyTypesSelectedThanChips(ColumnType type){
		return getCountOfType(type) == getChipCount();
	}
	
	/**
	 * Has user selected more columns for certain type than chips exists?
	 * 
	 * @param type
	 * @return
	 */
	public boolean isMoreTypesSelectedThanChips(ColumnType type){
		return getCountOfType(type) > getChipCount();
	}
	
	/**
	 * Does the column have a proper chip number set. A proper 
	 * chip number is number which is less or equal than total 
	 * chip count. If column type is <strong>annotation</strong> 
	 * the method return always true.
	 * 
	 * @param columnIndex
	 * @return
	 */
	public boolean isChipNumberSetProperly(int columnIndex){
		if(getColumnType(columnIndex).equals(ColumnType.UNUSED_LABEL)){
			return true;
		} else if(getColumnType(columnIndex).equals(ColumnType.ANNOTATION_LABEL)) { 
			return true;
		} else {
			int chip = getColumnChipNumber(columnIndex);
			return chip > 0 && chip <= getChipCount();
		}
	}
	
	/**
	 * Returns the count of correctly set columns for certain type. 
	 * Correctly set means that the type and chip number is correctly
	 * set
	 * 
	 * @param type
	 * @return
	 */
	public int getCountOfCorrectlySet(ColumnType type){
		int count = 0;
		for(int columnIndex = 0; columnIndex < columnCount; columnIndex++){
			ColumnType typeFromList = columns.get(columnIndex).getColumnType();
			if(typeFromList != null && typeFromList.equals(type)){
				if(isChipNumberSetProperly(columnIndex)){
					count++;
				}
			}
		}
		return count;
	}

	/**
	 * Sets new chip number to column. If chip number is already set to some other 
	 * column, removes the chip number from the column and sets number to the new 
	 * column.
	 * 
	 * @param columnIndex column index
	 * @param chipNumber new chip number
	 */
	public void setColumnChipNumber(int columnIndex, int chipNumber) {
		// Remove chip number if it is set before to some another column
		for(DataColumn column : getColumns()){
			if(column.getColumnType() == getColumnType(columnIndex) && 
					column.getChipNumber() == chipNumber &&
					column.getColumnIndex() != columnIndex){
				column.setChipNumber(-1);
				
				// Notify listeners
				fireChipNumberChangeEvent(new ChipNumberChangedEvent(this, Integer.MAX_VALUE, column.getColumnIndex()));
			}
		}
		
		// Set the chip number
		this.columns.get(columnIndex).setChipNumber(chipNumber);
		fireChipNumberChangeEvent(new ChipNumberChangedEvent(this, chipNumber, columnIndex));
	}

	/**
	 * Gets columns. The columns also includes the first column which 
	 * is the row number column.
	 * 
	 * @return columns and the row number column
	 */
	public List<DataColumn> getColumns() {
		return this.columns;
	}

	public void addColumnTypeChangeListener(ColumnTypeChangeListener l) {
		this.columnTypeListeners.add(l);
	}

	public void fireChipNumberChangeEvent(ChipNumberChangedEvent event) {
		for(ColumnTypeChangeListener listener : columnTypeListeners){
			listener.chipNumberChanged(event);
		}
		logger.debug("Fired event " + event.toString());
	}

	public void fireColumnTypeChangeEvent(ColumnTypeChangedEvent event) {
		for(ColumnTypeChangeListener listener : columnTypeListeners){
			listener.columnTypeChanged(event);
		}
		logger.debug("Fired event " + event.toString());
	}

	public void removeColumnTypeChangeListener(ColumnTypeChangeListener l) {
		this.columnTypeListeners.remove(l);
	}

	public void fireChipCountChangeEvent(ChipCountChangeEvent event) {
		for(ColumnTypeChangeListener listener : columnTypeListeners){
			listener.chipCountChanged(event);
		}
		logger.debug("Fired event " + event.toString());		
	}

	public void resetColumnTypes() {
		// Starts from 1 because of the row number column
		for(int columnIndex = 1; columnIndex < getColumnCount(); columnIndex++){
			setColumnType(columnIndex, ColumnType.UNUSED_LABEL, null);
		}
	}
	
}
