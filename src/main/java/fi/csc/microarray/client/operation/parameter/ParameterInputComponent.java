package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

/**
 * An abstract superclass for all the different input components specific
 * for different types of parameters.
 * 
 * @author Janne Käki
 *
 */
public abstract class ParameterInputComponent extends JPanel implements FocusListener {

	protected static final Color BG_VALID = Color.white;
	protected static final Color BG_INVALID = new Color(255, 190, 170);
	
	protected static final int PREFERED_WIDTH = 140;
	protected static final int PREFERED_HEIGHT = 23;
	protected static final Dimension PREFERRED_SIZE = new Dimension(
			ParameterInputComponent.PREFERED_WIDTH, 
			ParameterInputComponent.PREFERED_HEIGHT);
	
	protected static final int INPUT_IS_INITIALIZED = 0;
	protected static final int INPUT_IS_VALID = 1;
	protected static final int INPUT_IS_OUT_OF_BOUNDS = -1;
	protected static final int INPUT_IS_INCOMPREHENSIBLE = -2;
	
	private JLabel label = null;
	
	private ParameterPanel parentPanel;
	
	/**
	 * Initializes the features common to all parameter input components.
	 * To be used only by subclasses.
	 * 
	 * @param parent The parameter panel to which this component is to be placed.
	 */
	protected ParameterInputComponent(ParameterPanel parent) {
		super(new BorderLayout());
		this.parentPanel = parent;
	}
	
	/**
	 * @return The name label of this component.
	 */
	public JLabel getLabel() {
		if (label == null) {
			label = new JLabel(getParameter().getName());
		}
		return label;
	}
	
	public abstract JComponent getParameterComponent();
	
	/**
	 * @return The parent parameter panel to which this component is placed.
	 */
	public ParameterPanel getParentPanel() {
		return parentPanel;
	}
	
	/**
	 * @return The parameter that is edited via this component.
	 */
	public abstract Parameter getParameter();
	
	/**
	 * @return True if the user's input on this component is valid for the
	 * 		   annexed parameter, false if it's not.
	 */
	public abstract boolean inputIsValid();
	
	public void focusGained(FocusEvent e) {
	}

    public void focusLost(FocusEvent e){
    	getParentPanel().getOperationPanel().showOperationInfoText();
    }
}
