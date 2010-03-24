package fi.csc.microarray.databeans;

public class DataItemBase implements DataItem {

	protected String name;
	protected DataFolder parent;
	
	public String getName() {
		return name;
	}

	public void setName(String newName) {
		this.name = newName;
	}

	public DataFolder getParent() {
		return parent;
	}

	public void setParent(DataFolder newParent) {
		this.parent = newParent;
		
	}

}
