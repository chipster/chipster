package fi.csc.microarray.client.dataimport;

import java.util.ArrayList;
import java.util.List;

/**
 * Column pattern class for filling the columns with the same pattern from 
 * the beginning of the table to the end of the table.
 * 
 * @author mkoski, klemela
 *
 */
public class ColumnTypePattern {
	
	private List<ColumnType> pattern;
	private int patternStart;
	
	public ColumnTypePattern(List<ColumnType> pattern, int startIndex) {
		this.patternStart = startIndex;
		this.pattern = pattern;
	}
	
	public List<ColumnType> getPattern(){
		return pattern;
	}
	

	
	/** 
	 * @param i
	 * @return null if pattern doesn't cover asked index
	 */
	public ColumnType getColumnTypeForIndex(int i){
		if(i >= patternStart){
			return pattern.get((i - patternStart) % pattern.size());
		} else {
			return null;
		}
	}
	
	
	/**
	 * Tries to fiend pattern from the existing column types. Now this supports tree different use 
	 * cases: 
	 * 
	 * 1. 	There are marked column(s) repeating, possibly separated by constant amount of unused 
	 * 		columns. Also there can be some different types in the beginning, before the pattern
	 * 		starts. At least first column from the second repeat of the pattern has to be filled
	 * 
	 * 2. 	Pattern is just everything between first and last filled columns.
	 * 
	 * 3. 	There is no special pattern, just fill the rest of column with a one only type selected.
	 * 
	 * @param allColumns
	 * @return pattern
	 */
	public static ColumnTypePattern createColumnTypePatternFromAllColumns(List<DataColumn> allColumns){
		
		List<ColumnType> pattern = new ArrayList<ColumnType>();
		
		int patternEnd = -1;
		int start = -1;		
		
		//Start from end and find same column type as last filled
		for(int i = allColumns.size() - 1 ; i > 0; i--){ //Skip the first column (the row number column)
			ColumnType column = allColumns.get(i).getColumnType();
			if(column.equals(ColumnType.UNUSED_LABEL)){
				continue;
			} else if(patternEnd == -1){
				patternEnd = i;
				continue;
			} else if(column.equals(allColumns.get(patternEnd).getColumnType())){
				start = i + 1; //Don't repeat the found duplicate
				break;
			} 						
		}
		
		if(patternEnd == -1){//all columns are marked unused, give up
			return new ColumnTypePattern(pattern, Integer.MAX_VALUE);
		}
		
		//If no same type was found, search first any filled column
		if(start == -1){
			for(int i = 1; i < patternEnd; i++){
				ColumnType column = allColumns.get(i).getColumnType();
				if(column.equals(ColumnType.UNUSED_LABEL)){
					continue;
				} else {
					start = i;
					break;
				}
			}
		}
		
		//If still no start point found, just repeat the last type
		if(start == -1){
			start = patternEnd;
		}
				
		
		//Just collect the pattern
		for(int i = start; i <= patternEnd; i++){
			ColumnType column = allColumns.get(i).getColumnType();
			pattern.add(column);
		}
		
		return new ColumnTypePattern(pattern, start);
	}
	
	@Override
	public String toString(){
		StringBuffer pattern = new StringBuffer();
		for(ColumnType column : getPattern()){
			pattern.append(column.getIdentifier());
			pattern.append(" - ");
		}
		return "Column type pattern: " + pattern.toString();
	}
	
}
