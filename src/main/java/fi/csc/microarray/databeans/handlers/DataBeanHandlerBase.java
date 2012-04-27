package fi.csc.microarray.databeans.handlers;

import fi.csc.microarray.databeans.DataManager;

public abstract class DataBeanHandlerBase implements DataBeanHandler {

	protected DataManager dataManager;
	
	protected DataBeanHandlerBase(DataManager dataManager) {
		this.dataManager = dataManager;
	}
}
