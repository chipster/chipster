package fi.csc.microarray.client.operation;

import java.util.HashMap;
import java.util.LinkedList;

import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.databeans.DataBean;

public class OperationRecord {

	private String id;
	private String displayName;
	
	private String sourceCode;

	private LinkedList<Parameter> parameters;
	private HashMap<String, DataBean> inputs;
	
	
}
