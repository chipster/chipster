package fi.csc.microarray.databeans.handlers;

import java.util.HashSet;
import java.util.Set;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.StorageMethod;

public abstract class DataBeanHandlerBase implements DataBeanHandler {

	protected Set<StorageMethod> supportedDataBeanTypes;
	
	protected DataBeanHandlerBase(StorageMethod... supportedTypes) {
		this.supportedDataBeanTypes = new HashSet<StorageMethod>();
		for (StorageMethod type: supportedTypes) {
			this.supportedDataBeanTypes.add(type);
		}
	}

	protected void checkCompatibility(DataBean dataBean) {
		if (dataBean == null) {
			throw new IllegalArgumentException("DataBean is null.");
		}
		
		if (!supportedDataBeanTypes.contains(dataBean.getStorageMethod())) {
			throw new IllegalArgumentException("Unsupported DataBean type: " + dataBean.getStorageMethod());
		}
	}

}
