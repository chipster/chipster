package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.FocusEvent;
import java.text.DecimalFormat;

import javax.swing.AbstractSpinnerModel;
import javax.swing.JComponent;
import javax.swing.JFormattedTextField;
import javax.swing.JSpinner;
import javax.swing.JTextField;
import javax.swing.SpinnerModel;
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
		SpinnerModel model = new NullableSpinnerModel();
		model.setValue(param.getValue());
		this.spinner = new JSpinner(model);
		spinner.addFocusListener(this);
		spinner.setPreferredSize(ParameterInputComponent.PREFERRED_SIZE);	
		
		// The second parameter of NumberEditor constructor is number format
		// The string "0" means that it is a digit without any thousand separators
		// or decimals
		spinner.setEditor(new NullableSpinnerEditor(spinner, "0"));
		
		spinner.addChangeListener(this);
		field = ((JSpinner.DefaultEditor)spinner.getEditor()).getTextField();
		field.addCaretListener(this);
        field.setBackground(BG_VALID);
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
			NullableSpinnerModel numberModel =
				(NullableSpinnerModel) spinner.getModel();
			int value = numberModel.getNumber();
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
			if(spinner.getModel() instanceof NullableSpinnerModel){
				value = new Integer(((NullableSpinnerModel)spinner.getModel()).getNumber());
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

	/**
	 * A spinner model that can have an empty value. For example,
	 * when user has not entered anything. Normal SpinnerNumberModel
	 * would default to 0.
	 * 
	 * @author naktinis
	 */
	public class NullableSpinnerModel extends AbstractSpinnerModel {

	    protected String value = "";

	    public void setValue(Object o) {
	        if (o != null) {
	            value = o.toString();
	        } else {
	            value = "";
	        }
	        fireStateChanged();
	    }
	    
	    public Object getValue() {
	        return value;
	    }
	    
	    public Object getPreviousValue() {
	        Integer i = getNumber();
	        if (i == null) {
	            return "0";
	        } else { 
	            return "" + (i.intValue() - 1);
	        }
	    }
	    
	    public Object getNextValue() {
	        Integer i = getNumber();
	        if (i == null) {
	            return "0";
	        } else {
	            return "" + (i.intValue() + 1);
	        }
	    }

	    private Integer getNumber() {
	        try {
	            return new Integer(value);
	        } catch (NumberFormatException exc) {
	            return null;
	        }
	    }
	}
	
	/**
	 * Control for NullableSpinnerModel.
	 * 
	 * @author naktinis
	 */
	public class NullableSpinnerEditor extends JSpinner.DefaultEditor {
	    
	    private JSpinner spinner;
	    private DecimalFormat format;
	    private JTextField textField;

        public NullableSpinnerEditor(JSpinner spinner, String format) {
            super(spinner);
            
            this.spinner = spinner;
            this.format = new DecimalFormat(format);
            this.textField = getTextField();
            this.textField.setEditable(true);
            this.textField.setHorizontalAlignment(JTextField.RIGHT);
        }
        
        public DecimalFormat getFormat() {
            return format;
        }
        
        public NullableSpinnerModel getModel() {
            return (NullableSpinnerModel)spinner.getModel();
        }
	}

}