package fi.csc.microarray.databeans;

import fi.csc.microarray.databeans.Dataset.Link;

public class LinksChangedEvent extends DataChangeEvent {

	private Link type;
	private Dataset target;
	private boolean isCreation;

	public LinksChangedEvent(Dataset source, Dataset target, Link type, boolean isCreation) {
		super(source);
		this.target = target;
		this.type = type;
		this.isCreation = isCreation;
	}
	
	public Dataset getSource() {
		return (Dataset)getDataItem();
	}
	
	public Dataset getTarget() {
		return target;
	}
	
	public Link getType() {
		return type;
	}
	
	/**
	 * Returns true is this link was created, false if it was removed.
	 */
	public boolean isCreation() {
		return isCreation;
	}
}
