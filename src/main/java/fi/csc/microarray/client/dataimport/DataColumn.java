package fi.csc.microarray.client.dataimport;

/**
 * Class to store information of column's type and chip number.
 * 
 * @author mkoski
 *
 */
public class DataColumn {

	/**
	 * Column type
	 */
	private ColumnType columnType;
	
	/**
	 * Chip number. Negative value means that the chip number is not set
	 */
	private int chipNumber;

	private int columnIndex;
	
	private String originalName;
	
	public DataColumn(ColumnType columnType, int chipNumber, int columnIndex, String originalName) {
		this.columnType = columnType;
		this.chipNumber = chipNumber;
		this.columnIndex = columnIndex;
		this.originalName = originalName;
	}
	
	public int getChipNumber() {
		return chipNumber;
	}

	public void setChipNumber(int chipNumber) {
		this.chipNumber = chipNumber;
	}

	public ColumnType getColumnType() {
		return columnType;
	}
	
	public String getOriginalName(){
		return originalName;
	}

	public void setColumnType(ColumnType columnType) {
		if(columnType == null){
			columnType = ColumnType.UNUSED_LABEL;
		}
		this.columnType = columnType;
	}
	
	/**
	 * Creates an empty data column object and returns it. 
	 * Empty DataColumn is data column with type UNUSED_LABEL 
	 * and chip number -1
	 * 
	 * @return empty data column
	 */
	public static DataColumn createEmptyDataColumn(int columnIndex){
		return new DataColumn(ColumnType.UNUSED_LABEL, -1, columnIndex, null);
	}
	
	public static DataColumn createRowNumberColumn(int columnIndex){
		return new DataColumn(ColumnType.ROW_NUMBER, -1, columnIndex, null);
	}

	public int getColumnIndex() {
		return columnIndex;
	}

	public void setColumnIndex(int columnIndex) {
		this.columnIndex = columnIndex;
	}
	
	@Override
	public String toString(){
		return "DataColumn [Type: "+getColumnType().getIdentifier()+" Chip: " + getChipNumber() + " Index: " + getColumnIndex() + "]";
	}
	
}
