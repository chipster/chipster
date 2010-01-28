package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationDefinition.Suitability;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskEventListener;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.description.VVSADLSyntax.InputType;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * Operation is a concrete representative of an OperationDefinition, one with
 * unique parameters and a dataset associated to it. The parameters for a new
 * Operation are cloned from its definition object. Much other data is obtained
 * directly from the definition.
 * 
 * @author Janne KÃ¤ki, Aleksi Kallio
 * 
 */
public class Operation implements ExecutionItem {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(Operation.class);

	public static class DataBinding {

		private DataBean data;
		private String name;
		private InputType inputType;

		public DataBinding(DataBean data, String name, InputType inputType) {
			this.data = data;
			this.name = name;
			this.inputType = inputType;
		}

		public DataBean getData() {
			return data;
		}

		public String getName() {
			return name;
		}

		public InputType getInputType() {
			return inputType;
		}
	}

	private OperationDefinition definition;
	private LinkedList<DataBinding> bindings = new LinkedList<DataBinding>();
	private LinkedList<Parameter> parameters;
	private ResultListener resultListener;
	private DataFolder outputFolder;

	public static interface ResultListener {
		public void resultData(Iterable<DataBean> results);
		public void noResults();
	}

	/**
	 * Creates a new operation from the given definition and binds the given
	 * beans.
	 * 
	 * @param definition
	 * @param data
	 * @throws MicroarrayException
	 * @throws MicroarrayException
	 */
	public Operation(OperationDefinition definition, DataBean[] beans) throws MicroarrayException {
		logger.debug("created operation from definition " + definition.getName());
		this.definition = definition;
		this.parameters = Parameter.cloneParameters(definition.getParameters());
		bindInputs(beans);
	}

	/**
	 * Creates a clone of the given operation, with initially identical yet
	 * distinct attributes (such as parameters, which in turn are clones of the
	 * originals) but without a ParameterPanel.
	 * 
	 * @param o
	 *            The operation to be cloned.
	 * @throws MicroarrayException
	 */
	public Operation(Operation o) throws MicroarrayException {
		logger.debug("cloned operation from " + o.getName());
		this.definition = o.definition;
		this.bindings = o.bindings;
		this.parameters = Parameter.cloneParameters(o.getParameters());
	}

	/**
	 * @return The OperationDefinition from which this operation is derived.
	 */
	public OperationDefinition getDefinition() {
		return definition;
	}

	/**
	 * @return An array of this operation's parameters.
	 */
	public LinkedList<Parameter> getParameters() {
		return parameters;
	}

	/**
	 * @param beans operation input beans to bind
	 * @throws MicroarrayException 
	 */
	public void bindInputs(DataBean[] beans) throws MicroarrayException {
		this.bindings = definition.bindInputs(Arrays.asList(beans));

		// update bindings to parameters
		for (Parameter parameter : parameters) {
			if (bindings == null) {
				parameter.setDataBindings(new LinkedList<DataBinding>());
			} else {
				parameter.setDataBindings(this.bindings);
			}
		}
	}

	public void clearBindings() {
		this.bindings = null;
	}

	public List<DataBinding> getBindings() {
		return bindings;
	}

	/**
	 * @return The name of this operation (actually, of its definition).
	 */
	public String getName() {
		return definition.getName();
	}

	/**
	 * @return A String description of this operation's purpose (actually taken
	 *         directly from its definition).
	 */
	public String getDescription() {
		return definition.getDescription();
	}

	/**
	 * @return The name of the operation category to which this belongs.
	 */
	public String getCategoryName() {
		return definition.getCategory().getName();
	}

	/**
	 * @return The visualization color of the category
	 */
	public Color getCategoryColor() {
		return definition.getCategory().getColor();
	}

	/**
	 * @return A short String representation of this operation, consisting of
	 *         operation (definition's) name and the eventual dataset for which
	 *         it is assigned.
	 */
	public String toString() {
		return definition.getName();
	}

	/**
	 * Evaluates the suitability of this operation for the given dataset.
	 * 
	 * @param data
	 *            The dataset for which to evaluate.
	 * @return One of the OperationDefinition.Suitability enumeration, depending
	 *         on how suitable the operation is judged.
	 */
	public Suitability evaluateSuitabilityFor(Iterable<DataBean> data) {
		return definition.evaluateSuitabilityFor(data);
	}

	/**
	 * Starts the execution of this operation and creates a Job object for that
	 * matter. The job will be monitored by the listener given as parameter.
	 */
	public void execute(TaskEventListener listener) {
		TaskExecutor je = Session.getSession().getJobExecutor("client-job-executor");
		Task job = je.createTask(definition.getJobPhrase());
		job.addTaskEventListener(listener);

		try {
			for (DataBinding binding : bindings) {
				logger.debug("added input " + binding.name + " to job");
				job.addInput(binding.name, binding.data);
			}
			for (Parameter parameter : parameters) {
				logger.debug("operation has parameter value " + parameter.getValue());
				job.addParameter(parameter.getName(), parameter.getValue());
			}

			je.startExecuting(job);

		} catch (TaskException jex) {
			Session.getSession().getApplication().reportException(jex);
		}
	}

	/**
	 * Used in client scripting.
	 */
	public void setParameter(String name, Object value) {
		getParameter(name).setValue(value);
	}

	public void parseParameter(String name, String stringValue) {
		getParameter(name).parseValue(stringValue);
	}

	public Parameter getParameter(String name) {
		for (Parameter parameter : parameters) {
			if (parameter.getName().equals(name)) {
				return parameter;
			}
		}
		return null;
	}
	
	public int getColorCount() {
		return definition.getColorCount();
	}

	public ResultListener getResultListener() {
		return resultListener;
	}

	public void setResultListener(ResultListener resultListener) {
		this.resultListener = resultListener;
	}

	public DataFolder getOutputFolder() {
		return outputFolder;
	}

	public void setOutputFolder(DataFolder outputFolder) {
		this.outputFolder = outputFolder;
	}
}
