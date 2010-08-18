package fi.csc.microarray.client.operation;

import fi.csc.microarray.client.operation.OperationDefinition.Suitability;
import fi.csc.microarray.databeans.DataBean;

/**
 * A common interface for Operations, OperationDefinitions (which eventually
 * produce Operations) and Workflows (which consist of several Operations).
 * For now, this exists only to make things a bit smoother with the
 * operation selection lists.
 * 
 * @author Janne KÃ¤ki
 *
 */
public interface ExecutionItem {

	/**
	 * @return The name of this ExecutionItem.
	 */
	public String getID();
	
	public String getDisplayName();
	
	
	/**
	 * @return The category name of this ExecutionItem (either the name of an
	 * 		   OperationCategory or, for workflows, simply "Workflow".
	 */
	public String getCategoryName();
	
	/**
	 * @return A written description of this ExecutionItem's function
	 * 		   and purpose.
	 */
	public String getDescription();
	
	/**
	 * Evaluates the suitability of this ExecutionItem for the given dataset.
	 * 
	 * @param data The dataset for which to evaluate.
	 * @param currentSuitability Suitability known prior to execution
	 *        of this method. It can be overridden by this evaluation.
	 * @return One of the OperationDefinition.Suitability enumeration,
	 * 		   depending on how suitable the operation is judged.
	 */
	public Suitability evaluateSuitabilityFor(Iterable<DataBean> data,
	        Suitability currentSuitability);
}
