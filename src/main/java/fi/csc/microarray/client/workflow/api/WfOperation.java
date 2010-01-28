package fi.csc.microarray.client.workflow.api;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.Operation.ResultListener;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

public class WfOperation {

	private Operation operation;
	
	public WfOperation(WfOperationDefinition definition, WfDataBean[] datas) throws MicroarrayException {
		DataBean[] wrappedDatas = new DataBean[datas.length];
		for (int i = 0; i < datas.length; i++) {
			wrappedDatas[i] = datas[i].getWrapped();
		}
		this.operation = new Operation(definition.getWrapped(), wrappedDatas);
	}
	
	public void setParameter(String name, Object value) {
		this.operation.setParameter(name, value);
	}
	
	public void setResultListener(ResultListener resultListener) {
		this.operation.setResultListener(resultListener);
	}

	public Operation getWrapped() {
		return operation;
	}

}
