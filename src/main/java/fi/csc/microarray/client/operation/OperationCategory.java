package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.Vector;

/**
 * A class representing a category of similar Operations (or, to be more
 * accurate, OperationDefinitions). For example, there might be distinct
 * categories for all Preprocessing, Normalization, and Clustering operations.
 * Categories are used in the OperationChoicePanel to allow the selection of
 * operations by category.
 * 
 * @author Janne Kï¿½ki, Aleksi Kallio
 *
 */
public class OperationCategory {

	public static OperationCategory IMPORT_CATEGORY = new OperationCategory("Importing");
	public static OperationCategory USER_MODIFICATION_CATEGORY = new OperationCategory("Edited by user");
	
	/**
	 * Checks if the category is one of the predefined pseudo categories. Object identity is not required,
	 * as long as names match.
	 */
	public static boolean isPseudoCategory(OperationCategory category) {
		if (IMPORT_CATEGORY.getName().equals(category.getName())) {
			return true;
			
		} else if (USER_MODIFICATION_CATEGORY.getName().equals(category.getName())) {
			return true;
			
		} else {
			return false;
		}
	}
	private String name;
	private Vector<OperationDefinition> operations;
	private Color color;
	
	/**
	 * Creates a new, empty OperationCategory.
	 * 
	 * @param name Name of the category. Something general, such as
	 * 			   "Normalization" or "Clustering" for example - this is
	 * 			   currently used as the latter part of the name of any
	 * 			   operation in this category.
	 * @param type Type of the category (as the DataType of the datasets
	 * 			   produced by the operations of this category).
	 */
	public OperationCategory(String name) {
		this.name = name;
		this.operations = new Vector<OperationDefinition>();
	}
	
	/**
	 * @return The name of this category.
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * @return A Vector containing the operation definitions of this category.
	 */
	public Vector<OperationDefinition> getOperationList() {
		return operations;
	}
	
	/**
	 * Adds an operation to this operation category. The operation must be
	 * of same resulting DataType that this category has been labeled with.
	 * 
	 * @param operation The operation to be added.
	 * @return True if addition succeeded, false if it failed (in that case,
	 * 		   either the operation was null, of wrong resulting DataType,
	 * 		   or already present in this category.
	 */
	public void addOperation(OperationDefinition operation) {
		operations.add(operation);
	}
	
	
	/**
	 * Sets color for the category group which is used to visualize the workflow view.
	 * @param color
	 */
	public void setColor(Color c){
		color = c;
	}
	
	/**
	 * Gives the color for the category group which is used to visualize the workflow view.
	 * @return color
	 */
	public Color getColor(){
		return color;
	}
	
	/**
	 * @return Simply the name of this category.
	 */
	public String toString() {
		return name;
	}
}
