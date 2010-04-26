package fi.csc.microarray.databeans.handlers;

import java.util.HashSet;
import java.util.Set;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataBeanType;

public abstract class DataBeanHandlerBase implements DataBeanHandler {

	protected Set<DataBeanType> supportedDataBeanTypes;
	
	protected DataBeanHandlerBase(DataBeanType... supportedTypes) {
		this.supportedDataBeanTypes = new HashSet<DataBeanType>();
		for (DataBeanType type: supportedTypes) {
			this.supportedDataBeanTypes.add(type);
		}
	}

	protected void checkCompatibility(DataBean dataBean) {
		if (dataBean == null) {
			throw new IllegalArgumentException("DataBean is null.");
		}
		
		if (!supportedDataBeanTypes.contains(dataBean.getType())) {
			throw new IllegalArgumentException("Unsupported DataBean type: " + dataBean.getType());
		}
	}

}
