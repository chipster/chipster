package fi.csc.microarray.databeans;


/**
 * DataFolder is used to manage DataBean objects.
 * 
 * @see DataBean
 * @author Aleksi Kallio
 *
 */
public interface DataFolder extends DataItem {


	/**
	 * Add a bean or a subfolder to this folder.
	 */
	public void addChild(DataItem child);

	/**
	 * Remove a bean or a subfolder from this folder.
	 */
	public void removeChild(DataItem child);
	
	/**
	 * @return all beans and subfolders of this folder.
	 */
	public Iterable<DataItem> getChildren();
	
	/**
	 * @return the first subfolder with the given name
	 */
	public DataFolder getChildFolder(String name);
	
    /**
     * @return the count of beans and subfolders in this folder.
     */
    public int getChildCount();    
}