package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.LayoutManager;
import java.util.HashMap;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.exception.MicroarrayException;

/**
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
public abstract class ParameterPanel extends JPanel {

	protected Map<Parameter, ParameterInputComponent> paramMap;
	

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
			
		} else if (parameter instanceof SingleSelectionParameter) {
			return new SingleSelectionInputComponent((SingleSelectionParameter)parameter, this);
			
		} else if (parameter instanceof StringParameter) {
			return new StringInputComponent((StringParameter) parameter, this);
			
		} else {		
			throw new IllegalArgumentException("The given Parameter object, " + parameter.getName() + ", was not of recognized type!");
		}
	}
	
	public abstract void setMessage(String message, Color color);
}
