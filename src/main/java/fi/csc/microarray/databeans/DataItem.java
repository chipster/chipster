package fi.csc.microarray.databeans;

/**
 * A common superclass for main components of the databeans package:
 * DataBean and DataFolder.
 * 
 * @author Aleksi Kallio
 *
 */
public interface DataItem {
	
	/**
	 * @return superfolder of this item.
	 */
	public DataFolder getParent();
	
	public void setParent(DataFolder newParent);
	
	/**
	 * Name can be the name of a file/directory or some more user friendly name.
	 */
    public String getName();
    
    /**
     * @see #getName()
     */
    public void setName(String newName);

}