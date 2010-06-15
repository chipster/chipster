package fi.csc.microarray.client.workflow.api;

import fi.csc.microarray.databeans.Dataset;

public class WfDataBean {

	private Dataset dataBean;

	public WfDataBean(Dataset dataBean) {
		this.dataBean = dataBean;
	}

	public Dataset getWrapped() {
		return dataBean;
	}

}
