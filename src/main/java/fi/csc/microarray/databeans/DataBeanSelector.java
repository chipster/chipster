package fi.csc.microarray.databeans;

/**
 * DataBeanSelector is used to select which beans to traverse (search through) and
 * which beans to include in the result of the traversal.
 * 
 * @author Aleksi Kallio
 *
 */
public interface DataBeanSelector {
	
	/**
	 * Should this bean be included in the result of the traversal?
	 */
	public boolean shouldSelect(Dataset bean);
	
	/**
	 * Should this bean be traversed ie. links from it followed?
	 */
	public boolean shouldTraverse(Dataset bean);

}
