package fi.csc.microarray.client;

public class NameID {
	
	private String id;
	private String displayName;
	private String description;

	public NameID() {
		
	}

	public NameID(String id) {
		this(id, null, null);
	}
	
	public NameID(String id, String displayName, String description) {
		this.id = id;
		this.displayName = displayName;
		this.description = description;
		
	}


	public String getID() {
		return id;
	}




	public void setId(String id) {
		this.id = id;
	}



	public String getDisplayName() {
		return displayName;
	}




	public void setDisplayName(String displayName) {
		this.displayName = displayName;
	}




	public String getDescription() {
		return description;
	}




	public void setDescription(String description) {
		this.description = description;
	}




	public int hashCode() {
		return id.hashCode();
	}
}
