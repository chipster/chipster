package fi.csc.microarray.client.operation.parameter;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.parameter.SingleSelectionParameter.SelectionOption;
import fi.csc.microarray.description.SADLSyntax.ParameterType;
import fi.csc.microarray.exception.MicroarrayException;


/**
 * An abstract class representing all the different types of parameters
 * for the operations of this application.
 * 
 * @author Janne KÃ¤ki, akallio
 *
 */
public abstract class Parameter implements Cloneable {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(Parameter.class);

	/**
	 * Produces perfect clones of the given parameters - distinct entities
	 * with the same values and value limits. This method is mainly used
	 * when a new Operation is created out of an OperationDefinition.
	 * 
	 * @param parameters The array of parameters to be cloned.
	 * @return An array containing the clones.
	 */
	public static LinkedList<Parameter> cloneParameters(LinkedList<Parameter> parameters) {
		if (parameters == null) {
			return null;
		}
		
		LinkedList<Parameter> clones = new LinkedList<Parameter>();
		for (int i = 0; i < parameters.size(); i++) {
			clones.add((Parameter)parameters.get(i).clone());
		}
		return clones;
	}
	
	public void setDataBindings(List<DataBinding> bindings) throws MicroarrayException {
		// empty implementation, bindings are not usually needed
	}
	

	public static Parameter createInstance(String name, ParameterType type, String[] options, String description, String minValue, String maxValue, String initValue) {
		
		Parameter parameter = null;
		
		logger.debug("creating instance of parameter type " + type.name() + " called "+ name);
		
		switch (type) {
		case ENUM:
			int initIndex = (initValue != null ? Arrays.asList(options).indexOf(initValue) : 0);
			SelectionOption[] optionObjects = SingleSelectionParameter.SelectionOption.convertStrings(options);
			parameter = new SingleSelectionParameter(name, description, optionObjects, initIndex);
			break;
			
		case COLUMN_SEL:
			parameter = new ColnameParameter(name, description, initValue);
			break;

		case METACOLUMN_SEL:
			parameter = new MetaColnameParameter(name, description, initValue);
			break;

		case INPUT_SEL:
			parameter = new InputSelectParameter(name, description, initValue);
			break;

		case STRING:
			parameter = new StringParameter(name, description, initValue);
			break;
			
		case INTEGER:
		case DECIMAL:
		case PERCENT:
			
			// treat all numbers as double
			float min = (minValue == null ? Float.MIN_VALUE : Float.parseFloat(minValue));
			float max = (maxValue == null ? Float.MAX_VALUE : Float.parseFloat(maxValue));
			float init;
			if (initValue == null) {
				init = (min > 0F ? min : 0F); // default init is min or 0
			} else {
				init = Float.parseFloat(initValue);
			}
			
			switch (type) {
			case INTEGER:
				parameter = new IntegerParameter(name, description, (int)min, (int)max, (int)init);
				break;
				
			case DECIMAL:
				parameter = new DecimalParameter(name, description, min, max, init);
				break;
				
			case PERCENT:
				// put these to [0, 100]
				min = (min < 0F ? 0F : min);
				max = (max > 100F ? 100F : max);
				parameter = new PercentageParameter(name, description, (int)min, (int)max, (int)init);
				break;
			} 
			break;
			
		default:
			throw new IllegalArgumentException("unknown type " + type);
		}

		return parameter;
	}
	
	private String name;
	private String description;
	
	/**
	 * Constructor for initializing the name. Used by subclasses only.
	 * 
	 * @param name The name of this parameter.
	 */
	protected Parameter(String name, String description) {
		this.name = name;
		this.description = description;
	}
	
	/**
	 * @return The name of this parameter.
	 */
	public String getName() {
		return name;
	}
	
	/**
	 * A method used to access a parameter's value. Each subclass must
	 * implement this to return the appropriate type of parameter.
	 * 
	 * @return The value (be it a wrapped integer, double or something else)
	 * 		   of the parameter.
	 */
	public abstract Object getValue();
	
	/**
	 * Return the value in Java style formatting.
	 * 
	 * @return a representation of the value that can be inserted directly into Java code
	 */
	public abstract String getValueAsJava();
	
	
	/**
	 * A method used to set the value of a parameter. Each subclass must
	 * implement this and make it work appropriately.
	 * 
	 * @param newValue The new value of this parameter. Observe that each
	 * 				   subclass usually accepts only certain types of
	 * 				   parameters.
	 * @throws IllegalArgumentException If the given value cannot be accepted.
	 */
	public abstract void setValue(Object newValue) throws IllegalArgumentException;

	/**
	 * Tries to convert String representation to an Object that is valid value
	 * for this type of parameters and calls setValue to update the parameter
	 * values.
	 * 
	 * @see #setValue(Object)
	 * @param stringValue String representation of parameter value
	 * @throws IllegalArgumentException thrown if conversion fails
	 */
	public abstract void parseValue(String stringValue) throws IllegalArgumentException;

	/**
	 * Checks whether the given Object
	 * a) is of valid type to go as a value for this parameter, and
	 * b) is within the given value limits.
	 * If both tests succeed, the Object is judged valid, otherwise invalid.
	 * 
	 * @param valueObject The Object whose validity as a value is to be checked.
	 * @return True if argument can be accepted as a value for this parameter,
	 * 		   false if it cannot be.
	 */
	public abstract boolean checkValidityOf(Object valueObject);
	
	/**
	 * Clones this parameter by calling upon the superclass' (Object's)
	 * clone() method.
	 * 
	 * @return A clone (distinct entity with identical attributes) of
	 * 		   this parameter.
	 */
	public Object clone() {
		try {
			return super.clone();
		} catch (CloneNotSupportedException e) {
			return null;
		}
	}
	
	/**
	 * @return A String representation of this parameter. Subclasses shall
	 * 		   implement this and provide the appropriate texts.
	 */
	public abstract String toString();

	public String getDescription() {
		return description;
	}
}
