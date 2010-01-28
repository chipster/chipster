package fi.csc.microarray.client.workflow.api;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.wizard.ResultBlocker;
import fi.csc.microarray.databeans.DataBean;

public class WfResultBlocker extends ResultBlocker {

	public WfResultBlocker() {
		super();
	}

	public WfResultBlocker(int enforcedResultCount) {
		super(enforcedResultCount);
	}

	public List<WfDataBean> getWorkflowResults() {
		List<DataBean> results = getResults();
		LinkedList<WfDataBean> wfResults = new LinkedList<WfDataBean>();

		for (DataBean result : results) {
			wfResults.add(new WfDataBean(result));
		}
		
		return wfResults;
	}
}
