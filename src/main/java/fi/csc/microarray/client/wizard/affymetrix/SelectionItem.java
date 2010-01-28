package fi.csc.microarray.client.wizard.affymetrix;

public class SelectionItem {
	private String publicValue = null;
	private String internalValue = null;
	
	public SelectionItem(String publicVal, String internalVal) {
		publicValue = publicVal;
		internalValue = internalVal;
	}
	
	public String getPublicValue() {
		return publicValue;
	}
	
	public String getInternalValue() {
		return internalValue;
	}
}
