package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.FocusEvent;

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JSlider;
import javax.swing.SwingConstants;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * A component for controlling the value of a PercentageParameter.
 * This is done with a JSlider, accompanied with a label displaying the
 * percentage in numbers. For now, the value can be adjusted with the
 * accuracy of 1 (although for some reason, 99% seems inaccessible..).
 * Sliding the value out of the given limits (which can be something
 * other than 0 and 100 percent) is technically allowed, but the label
 * will display the value in bright red.
 * 
 * @author Janne KÃ¤ki
 *
 */
@SuppressWarnings("serial")
public class PercentageInputComponent extends ParameterInputComponent
								      implements ChangeListener {

	private JSlider slider;
	private JLabel numberLabel;
	
	private PercentageParameter param;
	private int state;
	
	/**
	 * Creates a new PercentageInputComponent.
	 * 
	 * @param param The PercentageParameter to be controlled.
	 * @param enabled 
	 * @param parameterPanel The ParameterPanel to which this component is to
	 * 				 be placed.
	 */
	public PercentageInputComponent(
			PercentageParameter param, ParameterPanel parameterPanel) {
		super(parameterPanel);
		this.param = param;
		this.state = ParameterInputComponent.INPUT_IS_INITIALIZED;
		/**
		 * The eventual minimum and maximum values of the percentage parameter
		 * may be something between 0 and 100, but the slider will always show
		 * a scale from 0 to 100, even if some values would be illegal.
		 * The illegality is (hopefully) clearly indicated.
		 */
		int initValue = 0; 
		if (param.getIntegerValue() != null) {
		    initValue = param.getIntegerValue().intValue();
		}
		this.slider = new JSlider(JSlider.HORIZONTAL, 0, 100, initValue);
		slider.setMajorTickSpacing(50);
		slider.setMinorTickSpacing(10);
		slider.setPaintTicks(true);
		slider.setPaintLabels(false);
		slider.addChangeListener(this);
		slider.addFocusListener(this);
		// Normal width minus label width
		slider.setPreferredSize(new Dimension(ParameterInputComponent.PREFERED_WIDTH - 40, 40));
		
		this.numberLabel = new JLabel(
				param.getIntegerValue() + "%", SwingConstants.CENTER);
		numberLabel.setPreferredSize(new Dimension(40, 30));
		numberLabel.setBorder(BorderFactory.createEmptyBorder(0, 10, 0, 0));
		
		this.add(slider, BorderLayout.WEST);
		this.add(numberLabel, BorderLayout.EAST);
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
		if (e.getSource() == slider) {
			Integer value = slider.getValue();
			numberLabel.setText(value + "%");
			// slider.setToolTipText(value + "%");
			if (param.checkValidityOf(value) == true) {
				param.setValue(value);
				setState(ParameterInputComponent.INPUT_IS_VALID);
			} else {
				setState(ParameterInputComponent.INPUT_IS_OUT_OF_BOUNDS);
			}
		}
	}
	
	private void setState(int newState) {
		// if (this.state != newState) {
		String message = null;
		this.state = newState;
		switch (state) {
		case ParameterInputComponent.INPUT_IS_VALID:
			numberLabel.setForeground(Color.black);
			message = param.getDescription();
			getParentPanel().setMessage(message, Color.black);
			break;
		case ParameterInputComponent.INPUT_IS_OUT_OF_BOUNDS:
			numberLabel.setForeground(Color.red);
			message =
				"Value for " + param.getName() + " must be between " +
				param.getMinValue() + "% and " + param.getMaxValue() + "%.";
			getParentPanel().setMessage(message, Color.red);
		}
	}

	@Override
	public JComponent getParameterComponent() {
		return slider;
	}

	public void focusGained(FocusEvent e) {
		getParentPanel().setMessage(param.getDescription(), Color.black);
	}


}