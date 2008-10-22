package fi.csc.microarray.client.dataimport.trimmer;

import java.util.List;
import java.util.Stack;
import java.util.Vector;


/**
 * Class for import screen find-replace and undo functionality. The data 
 * modification operations are stored in to a stack and they can be 
 * added and removed (undo) easily. The flaging functionality is also 
 * made with this class.
 * 
 * @author mkoski
 *
 */
public class DataTrimmer {

	private Stack<DataTrimmingOperation> operations;
	private Vector<Integer> ignoreColumns;
	
	public DataTrimmer(){
		operations = new Stack<DataTrimmingOperation>();
		ignoreColumns = new Vector<Integer>();
	}
	
	public void pushOperation(DataTrimmingOperation operation){
		operations.push(operation);
	}
	
	public void pushOperations(List<DataTrimmingOperation> operationList){
		operations.addAll(operationList);
	}
	
	public void popOperation(){
		operations.pop();
	}
	
	public int getOperationCount(){
		return operations.size();
	}
	
	public void clear(){
		operations.clear();
	}
	
	/**
	 * Column which will be ignored
	 *
	 */
	public void addIgnoreColumnNumber(int columnIndex){
		ignoreColumns.add(columnIndex);
	}
	
	/**
	 * Does the data modification with operations currently on the stack
	 * 
	 * @param stringToTrim string to be trimmed
	 * @param column column number that the string represents
	 * @return modified string
	 */
	public String doTrimming(String stringToTrim, int column){
		String stringToReturn = stringToTrim;
		for(DataTrimmingOperation operation : operations){
			if(!ignoreColumns.contains(column) && ((column == operation.getColumnIndex() || (operation.getColumnIndex() == DataTrimmingOperation.ALL_COLUMNS)))){
				stringToReturn = operation.doTrimming(stringToReturn);
			}
		}
		return stringToReturn;
	}
}
