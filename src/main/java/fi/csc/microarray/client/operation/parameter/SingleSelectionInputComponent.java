package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;

import javax.swing.JComponent;

import org.apache.log4j.Logger;

/**
 * A component for controlling the value of a EnumParameter.
 * This is done by using a JComboBox which displays all the possible
 * value objects. Thus, no erroneous selection is possible as the user
 * cannot give any undefined input (although a certain value of some other
 * parameter might make a choice semantically impossible - maybe this
 * possibility should be also taken into account).
 * 
 * @author Janne Käki
 *
 */
@SuppressWarnings("serial")
public class SingleSelectionInputComponent extends ParameterInputComponent
										   implements ActionListener {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(SingleSelectionInputComponent.class);

	private SteppedComboBox choiceBox;
	
	private EnumParameter param;
	
	/**
	 * Creates a new SingleSelectionInputComponent.
	 * 
	 * @param param The EnumParameter to be controlled.
	 * @param enabled 
	 * @param parameterPanel The ParameterPanel to which this component is to
	 * 				 be placed.
	 */
	public SingleSelectionInputComponent(
			EnumParameter param, ParameterPanel parameterPanel) {
		super(parameterPanel);
		this.param = param;
		this.choiceBox = new SteppedComboBox(param.getOptions());
		Dimension preferredSize = choiceBox.getPreferredSize();
		choiceBox.setPreferredSize(ParameterInputComponent.PREFERRED_SIZE);
		choiceBox.setPopupWidth(preferredSize.width);
		choiceBox.setBackground(Color.white);
		
		if (param.getSelectedOptions().size() > 0) {
		    // This is a single selection list so there's only one default value
		    choiceBox.setSelectedItem(param.getSelectedOptions().get(0));
		}
        
		choiceBox.addActionListener(this);
		choiceBox.addFocusListener(this);
		this.add(choiceBox, BorderLayout.CENTER);
	}
	
	@Override
	public Parameter getParameter() {
		return param;
	}
	
	@Override
	public boolean inputIsValid() {
		/**
		 * Currently we assume that any of the given options is
		 * automatically valid as an input. At least here the user
		 * has not the chance to input whatever arbitrary data they
		 * wish, as is the case with text fields.
		 */
		return true;
	}
	
	public void actionPerformed(ActionEvent e) {
		/**
		 * The indexing of the JComboBox is automagically the same as
		 * in the value item array of the parameter.
		 */
		param.setValue(choiceBox.getSelectedItem());
		logger.debug("set selection to " + choiceBox.getSelectedItem());
		getParentPanel().setMessage(param.getDescription(), Color.black);
	}

	public void focusGained(FocusEvent e) {
		getParentPanel().setMessage(param.getDescription(), Color.black);
	}

	@Override
	public JComponent getParameterComponent() {
		return choiceBox;
	}
}