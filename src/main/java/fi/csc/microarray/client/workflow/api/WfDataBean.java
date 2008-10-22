package fi.csc.microarray.client.workflow.api;

import fi.csc.microarray.databeans.DataBean;

public class WfDataBean {

	private DataBean dataBean;

	public WfDataBean(DataBean dataBean) {
		this.dataBean = dataBean;
	}

	public DataBean getWrapped() {
		return dataBean;
	}

}
