package fi.csc.microarray.client.operation.parameter;

import java.awt.Color;
import java.awt.LayoutManager;

import javax.swing.JPanel;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * The class ParameterPanel class hierarchy is a mess.
 * 
 * ParameterPanel is the GUI component for showing the parameters for an
 * operation, along with a contextual help textbox showing information about
 * the currently selected parameter and eventual errors in user input. 
 *  A new ParameterPanel will be created every time it is
 * to be shown, whether the corresponding Operation (derived from an
 * OperationDefinition) already exists or not. At least for now. Thought
 * it would save some memory.
 * 
 * @author Janne KÃ¤ki, Aleksi Kallio, Petri KlemelÃ¤
 *
 */
@SuppressWarnings("serial")
public abstract class ParameterPanel extends JPanel {

	public ParameterPanel(Operation operation) {
	}
	
	public ParameterPanel(Operation operation, LayoutManager layoutManager) {
		super(layoutManager);
	}
	
	/**
	 * Creates the appropriate type of parameter input component, depending
	 * of the type of the given parameter Parameter (pun intended).
	 * 
	 * @param parameter The parameter for which an input component is to be made.
	 * @return A GUI component that can be used by user to set the value
	 * 		   of the given parameter.
	 * @throws MicroarrayException 
	 */
	protected ParameterInputComponent createInputComponent(Parameter parameter) throws MicroarrayException {
		if (parameter instanceof IntegerParameter) {
			return new IntegerInputComponent((IntegerParameter) parameter, this);
			
		} else if (parameter instanceof DecimalParameter) {
			return new DecimalInputComponent((DecimalParameter) parameter, this);
			
		} else if (parameter instanceof PercentageParameter) {
			return new PercentageInputComponent((PercentageParameter) parameter, this);
			
		} else if (parameter instanceof DataSelectionParameter) {
			return new SingleSelectionInputComponent((DataSelectionParameter)parameter, this);
			
		} else if (parameter instanceof EnumParameter) {
		    EnumParameter enumParam = (EnumParameter) parameter;
		    if (enumParam.getMaxCount() > 1) {
		        // List with multiple selections
		        return new MultipleSelectionInputComponent(enumParam, this);
		    } else {
		        // List with single selection
		        return new SingleSelectionInputComponent(enumParam, this);
		    }
		} else if (parameter instanceof StringParameter) {
			return new StringInputComponent((StringParameter) parameter, this);
			
		} else {		
			throw new IllegalArgumentException("The given Parameter object, " + parameter.getID() + ", was not of recognized type!");
		}
	}
	
	public abstract void setMessage(String message, Color color);
}
