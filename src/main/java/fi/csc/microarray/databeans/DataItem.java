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
	public abstract DataFolder getParent();
	
	/**
	 * Name can be the name of a file/directory or some more user friendly name.
	 */
    public abstract String getName();
    
    /**
     * @see #getName()
     */
    public abstract void setName(String newName);
    
    public abstract void setParent(DataFolder newParent);

}