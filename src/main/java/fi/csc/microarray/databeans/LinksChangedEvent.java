package fi.csc.microarray.databeans;

import fi.csc.microarray.databeans.DataBean.Link;

public class LinksChangedEvent extends DataChangeEvent {

	private Link type;
	private DataBean target;
	private boolean isCreation;

	public LinksChangedEvent(DataBean source, DataBean target, Link type, boolean isCreation) {
		super(source);
		this.target = target;
		this.type = type;
		this.isCreation = isCreation;
	}
	
	public DataBean getSource() {
		return (DataBean)getDataItem();
	}
	
	public DataBean getTarget() {
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
