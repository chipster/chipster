package fi.csc.microarray.client.workflow.api;

import fi.csc.microarray.client.selection.DataSelectionManager;

public class WfSelectionManager {

	private DataSelectionManager selectionManager;

	public WfSelectionManager(DataSelectionManager selectionManager) {
		this.selectionManager = selectionManager; 
	}

	public WfDataBean getSelectedDataBean() {
		return new WfDataBean(selectionManager.getSelectedDataBean());
	}

}
