package fi.csc.microarray.client.workflow.api;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.tasks.ResultBlocker;
import fi.csc.microarray.databeans.Dataset;

public class WfResultBlocker extends ResultBlocker {

	public WfResultBlocker() {
		super();
	}

	public WfResultBlocker(int enforcedResultCount) {
		super(enforcedResultCount);
	}

	public List<WfDataBean> getWorkflowResults() {
		List<Dataset> results = getResults();
		LinkedList<WfDataBean> wfResults = new LinkedList<WfDataBean>();

		for (Dataset result : results) {
			wfResults.add(new WfDataBean(result));
		}
		
		return wfResults;
	}
}
