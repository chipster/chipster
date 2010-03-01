package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.FocusEvent;

import javax.swing.JComponent;
import javax.swing.JFormattedTextField;
import javax.swing.JSpinner;
import javax.swing.SpinnerModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * A component for controlling the value of an IntegerParameter.
 * This is done with a text field accompanied by a JSpinner (component with
 * arrows up and down). It will allow inserting of (and spinning to) illegal
 * values, but will indicate this with a red background. If the value cannot
 * be parsed to an integer at all, the input text will also be bright red.
 * 
 * @author Janne KÃ¤ki
 *
 */
@SuppressWarnings("serial")
public class IntegerInputComponent extends ParameterInputComponent
								   implements CaretListener, ChangeListener {

	private JSpinner spinner;
	private JFormattedTextField field;
	
	private IntegerParameter param;
	private int state;
	
	/**
	 * Creates a new IntegerInputComponent.
	 * 
	 * @param param The IntegerParameter to be controlled.
	 * @param parameterPanel The ParameterPanel to which this component is to
	 * 				 be placed.
	 */
	public IntegerInputComponent(
			IntegerParameter param, ParameterPanel parameterPanel) {
		super(parameterPanel);
		this.param = param;
		this.state = ParameterInputComponent.INPUT_IS_INITIALIZED;
		SpinnerModel model = new SpinnerNumberModel();
		model.setValue(param.getIntegerValue());
		this.spinner = new JSpinner(model);
		spinner.addFocusListener(this);
		spinner.setPreferredSize(ParameterInputComponent.PREFERRED_SIZE);	
		
		// The second parameter of NumberEditor constructor is number format
		// The string "0" means that it is a digit without any thousand separators
		// or decimals
		spinner.setEditor(new JSpinner.NumberEditor(spinner, "0"));
		
		spinner.addChangeListener(this);
		field = ((JSpinner.DefaultEditor)spinner.getEditor()).getTextField();
		field.addCaretListener(this);
		this.add(spinner, BorderLayout.CENTER);
	}	
	
	@Override
	public Parameter getParameter() {
		return param;
	}
	
	@Override
	public boolean inputIsValid() {
		if (state == ParameterInputComponent.INPUT_IS_INITIALIZED ||
				state == ParameterInputComponent.INPUT_IS_VALID) {
			return true;
		} else {
			return false;
		}
	}
	
	public void stateChanged(ChangeEvent e) {
		if (e.getSource() == spinner) {
			SpinnerNumberModel numberModel =
				(SpinnerNumberModel) spinner.getModel();
			int value = numberModel.getNumber().intValue();
			if (param.checkValidityOf(value) == true) {
				param.setValue(value);
				setState(ParameterInputComponent.INPUT_IS_VALID);
			} else {
				setState(ParameterInputComponent.INPUT_IS_OUT_OF_BOUNDS);
			}
		}
	}
	
	public void caretUpdate(CaretEvent e) {
		try {
			Integer value = null;
			if(spinner.getModel() instanceof SpinnerNumberModel){
				value = new Integer(((SpinnerNumberModel)spinner.getModel()).getNumber().intValue());
			}
			if (param.checkValidityOf(value) == true) {
				setState(ParameterInputComponent.INPUT_IS_VALID);
				param.setValue(value);
			} else {
				setState(ParameterInputComponent.INPUT_IS_OUT_OF_BOUNDS);
			}
		} catch (NumberFormatException nfe) {
			setState(ParameterInputComponent.INPUT_IS_INCOMPREHENSIBLE);
		}
	}
	
	private void setState(int newState) {
		String message = null;
		this.state = newState;
		switch (state) {
		case ParameterInputComponent.INPUT_IS_VALID:
			field.setBackground(ParameterInputComponent.BG_VALID);
			field.setForeground(Color.black);
			message = param.getDescription();
			getParentPanel().setMessage(message, Color.black);
			break;
		case ParameterInputComponent.INPUT_IS_OUT_OF_BOUNDS:
			field.setBackground(ParameterInputComponent.BG_INVALID);
			field.setForeground(Color.black);
			message =
				"Value for " + param.getName() + " must be between " +
				param.getMinValue() + " and " + param.getMaxValue() + ".";
			getParentPanel().setMessage(message, Color.red);
			break;
		case ParameterInputComponent.INPUT_IS_INCOMPREHENSIBLE:
			field.setBackground(ParameterInputComponent.BG_INVALID);
			field.setForeground(Color.red);
			message =
				"Value for " + param.getName() + " must be a valid integer.";
			getParentPanel().setMessage(message, Color.red);
		}
	}

	@Override
	public JComponent getParameterComponent() {
		return spinner;
	}

	public void focusGained(FocusEvent e) {
		getParentPanel().setMessage(param.getDescription(), Color.black);
	}



}