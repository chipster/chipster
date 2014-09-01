package fi.csc.microarray.client.workflow;

import java.io.IOException;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.InputRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.operation.parameter.DataSelectionParameter;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;

public class WorkflowWriter {

	public static String generateVersionHeaderLine() {
		return "// VERSION " + WorkflowManager.WORKFLOW_VERSION + " (do not remove this)\n";
	}

	private LinkedList<String> writeWarnings = new LinkedList<String>();
	private HashMap<DataBean, String> resultIdMap;
	private boolean used = false;

	/**
	 * Saves workflow from currently selected databean. Method is synchronised
	 * because it uses state of the manager during the generation process.
	 */
	public synchronized StringBuffer writeWorkflow(DataBean root) throws IOException {

		if (this.used) {
			throw new IllegalStateException("writer cannot be reused");
		}

		// initialise recursion
		this.resultIdMap = new HashMap<DataBean, String>();

		// generate header
		StringBuffer script = new StringBuffer("");
		generateHeader(script, root);

		// do recursion
		generateRecursively(script, root);

		// mark this writer dirty
		this.used = true;

		return script;
	}

	/**
	 * FIXME session refactoring, start using OperationRecords.
	 * @param script
	 * @param bean
	 */
	private void generateRecursively(StringBuffer script, DataBean bean) {

		// dig all operations that were used to produce children of this bean
		HashSet<OperationRecord> operationRecords = new HashSet<OperationRecord>();
		for (DataBean derived : bean.getLinkSources(Link.DERIVATION)) {
			operationRecords.add(derived.getOperationRecord());
		}

		// generate operations leading to derived beans
		for (OperationRecord operationRecord : operationRecords) {
			LinkedList<DataBean> results = new LinkedList<DataBean>();
			// collect datas that were produced by this operation
			for (DataBean derived : bean.getLinkSources(Link.DERIVATION)) {
				if (derived.getOperationRecord() == operationRecord) {
					results.add(derived);
				}
			}

			generateStep(script, operationRecord, results);
		}

		// continue recursion from derived beans
		for (DataBean derived : bean.getLinkSources(Link.DERIVATION)) {
			generateRecursively(script, derived);
		}
	}

	/**
	 * 
	 * @param results must be in the same order as they were produced by the operation!
	 */
	private void generateStep(StringBuffer script, OperationRecord operationRecord, LinkedList<DataBean> results) {
		
		OperationDefinition toolDefinition = Session.getSession().getApplication().getOperationDefinitionBestMatch(operationRecord.getNameID().getID(), operationRecord.getModule(), operationRecord.getCategoryName());
		if (toolDefinition == null) {
			// FIXME writeWarnings
		}
		
		StringBuffer dataString = new StringBuffer("\ndatas = new WfDataBean[] {\n");
		boolean first = true;
		for (InputRecord inputRecord : operationRecord.getInputRecords()) {
			
			// skip phenodata as it is bound automatically
			if (inputRecord.getNameID().getID().equals("phenodata.tsv")) {
				continue; 
			}
			if (!first) {
				dataString.append(",\n");
			} else {
				first = false;
			}

			
			String name = resultIdMap.get(inputRecord.getValue());
			if (name == null) {
				// we cannot handle this, too complicated structure
				writeWarnings.add("Tool " + operationRecord.getFullName() + " was skipped because it combines multiple workflow branches.");
				return; // skip this branch, nothing was written to script yet
			}
			dataString.append("  " + name);
		}
		dataString.append("\n};\n");
		script.append(dataString);

		script.append("op = new WfOperation(app.getOperationDefinitionBestMatch(\"" + operationRecord.getNameID().getID() + "\", \"" + operationRecord.getModule() + "\", \"" + operationRecord.getCategoryName() + "\"), datas);\n");

		for (ParameterRecord parameterRecord : operationRecord.getParameters()) {
			// Write code that sets the value only when value is not empty
			if (parameterRecord.getValue() != null && !parameterRecord.getValue().equals("")) {	
				Parameter parameter = toolDefinition.getParameter(parameterRecord.getNameID().getID());
				if (parameter != null) {
					parameter = (Parameter) parameter.clone();
					if (parameter instanceof DataSelectionParameter) {
						((DataSelectionParameter)parameter).parseValueAndSetWithoutChecks(parameterRecord.getValue());
					} else {
						parameter.parseValue(parameterRecord.getValue());
					}
					script.append("op.setParameter(\"" + parameter.getID() + "\", " + parameter.getValueAsJava() + ");\n");
				}
			}
		}

		script.append("opBlocker = new WfResultBlocker();\n");
		script.append("op.setResultListener(opBlocker);\n");
		script.append("app.executeOperation(op);\n");
		script.append("opBlocker.blockUntilDone();\n");

		int i = -1;
		for (DataBean result : results) {
			i++;
			// here we demand that the order of results is correct!
			String name = "data" + this.resultIdMap.size();
			this.resultIdMap.put(result, name);
			script.append(name + " = opBlocker.getWorkflowResults().get(" + i + ");\n");
		}
	}

	private void generateHeader(StringBuffer script, DataBean root) {
		// write info and imports
		script.append(generateVersionHeaderLine());
		script.append("/* \n" + "  BeanShell workflow script for " + ApplicationConstants.TITLE + "\n" + "  Generated by " + System.getProperty("user.name") + " at " + new Date().toString() + "\n" + "*/\n");
		script.append("\n");
		script.append("import fi.csc.microarray.client.workflow.api.*;\n");
		script.append("\n");
		
		// write datas
		script.append("data0 = app.getSelectionManager().getSelectedDataBean();\n");
		this.resultIdMap.put(root, "data0");
	}

	/**
	 * List operations that we had to skip.
	 */
	public List<String> writeWarnings() {
		return writeWarnings ;
	}
}
