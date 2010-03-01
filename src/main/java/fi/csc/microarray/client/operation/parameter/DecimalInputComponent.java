package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.FocusEvent;

import javax.swing.JComponent;
import javax.swing.JTextField;
import javax.swing.event.CaretEvent;
import javax.swing.event.CaretListener;

/**
 * A component for controlling the value of a DecimalParameter.
 * This is done with a simple text field. It will allow inserting of illegal
 * values, but will indicate this with a red background. If the value cannot
 * be interpreted as a decimal number at all (cannot be parsed to a double),
 * the input text will also be bright red.
 * 
 * @author Janne KÃ¤ki
 *
 */
@SuppressWarnings("serial")
public class DecimalInputComponent extends ParameterInputComponent
								   implements CaretListener {

	private JTextField field;
	
	private DecimalParameter param;
	private int state;
	
	/**
	 * Creates a new DecimalInputComponent.
	 * 
	 * @param param The DecimalParameter to be controlled.
	 * @param enabled 
	 * @param parameterPanel The ParameterPanel to which this component is to
	 * 				 be placed.
	 */
	public DecimalInputComponent(
			DecimalParameter param, ParameterPanel parameterPanel) {
		super(parameterPanel);
		this.param = param;
		this.field = new JTextField();
		field.setPreferredSize(ParameterInputComponent.PREFERRED_SIZE);
		field.setText("" + param.getDecimalValue());
		field.addCaretListener(this);
		field.addFocusListener(this);
		this.add(field, BorderLayout.CENTER);
		this.state = ParameterInputComponent.INPUT_IS_INITIALIZED;
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
	
	public void caretUpdate(CaretEvent e) {
		try {
			double value = Double.parseDouble(field.getText());
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
		// if (this.state != newState) {
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
				"Value for " + param.getName() + " must be a valid " +
				"decimal number. Use a point as the decimal separator.";
			getParentPanel().setMessage(message, Color.red);
		}
	}

	@Override
	public JComponent getParameterComponent() {
		return field;
	}
	
	public void focusGained(FocusEvent e) {
		getParentPanel().setMessage(param.getDescription(), Color.black);
	}

	
}