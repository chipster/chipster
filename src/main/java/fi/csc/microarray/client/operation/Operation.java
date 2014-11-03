package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.OperationDefinition.Suitability;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * Operation is a concrete representative of an OperationDefinition, one with
 * unique parameters and a dataset associated to it. The parameters for a new
 * Operation are cloned from its definition object. Much other data is obtained
 * directly from the definition.
 * 
 * @author Janne Käki, Aleksi Kallio, hupponen
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
		logger.debug("created operation from definition " + definition.getID());
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
		logger.debug("cloned operation from " + o.getID());
		this.definition = o.definition;
		if (o.bindings != null) {
			this.bindings = o.bindings;
		} else {
			this.bindings = new LinkedList<DataBinding>(); 
		}
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
		this.bindings = definition.bindInputs(Arrays.asList(beans)).bindings;
		if (bindings == null) {
			bindings = new LinkedList<DataBinding>();
		}
		
		// update bindings to parameters
		for (Parameter parameter : parameters) {
			parameter.setDataBindings(this.bindings);
		}
	}

	public void clearBindings() {
		this.bindings.clear();
	}

	public List<DataBinding> getBindings() {
		return bindings;
	}
	
	/**
	 * 
	 * @param name
	 * @return null if a binding with the given name does not exist
	 */
	public DataBinding getBinding(String name) {
		for (DataBinding binding : bindings) {
			if (binding.getName().equals(name)) {
				return binding;
			}
		}
		return null;
	}
	
	/**
	 * Set new input file bindings.
	 */
	public void setBindings(List<DataBinding> currentBindings) {
	    if (currentBindings != null) {
	    	this.bindings = new LinkedList<>(currentBindings);
	    } else {
	    	this.bindings = new LinkedList<DataBinding>();	
	    }
	}

	/**
	 * @return The name of this operation (actually, of its definition).
	 */
	public String getID() {
		return definition.getID();
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
		return definition.getID();
	}

	/**
	 * Evaluates the suitability of this operation for the given dataset
	 * and current parameter values.
	 * 
	 * @param data
	 *            The dataset for which to evaluate.
	 * @return One of the OperationDefinition.Suitability enumeration, depending
	 *         on how suitable the operation is judged.
	 */
	@Override
	public Suitability evaluateSuitabilityFor(Iterable<DataBean> data, List<DataBinding> bindings) {
	    
	    // Check parameters and inputs
	    Suitability parameterSuitability = 
	        OperationDefinition.evaluateParameterSuitability(getParameters());
	    Suitability inputSuitability = definition.evaluateSuitabilityFor(data, bindings);
	    
	    // Return ok suitability only when both are ok
	    if (inputSuitability.isOk()) {
	    	return parameterSuitability;
	    }  else {
	    	return inputSuitability;
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
			if (parameter.getID().equals(name)) {
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
	
	public String getDisplayName() {
		return definition.getDisplayName();
	}
}
