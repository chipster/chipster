package fi.csc.microarray.client.operation.parameter;

import java.awt.Color;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationPanel;
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
public class ImportParameterPanel extends ParameterPanel {

	private OperationPanel parent;
	
	/**
	 * Creates a new ParameterPanel for the given operation.
	 * 
	 * @param operation The operation which is to be started from this panel.
	 * @param client The client that is to be informed when a job is done.
	 * @throws MicroarrayException 
	 */
	public ImportParameterPanel(Operation operation, OperationPanel parent) throws MicroarrayException {
		super(operation);
		this.parent = parent;		
	}
	
	
	/**
	 * Sets the message of the (bottom left) infobox of this parameter panel.
	 * 
	 * @param message The text to be set.
	 * @param color The font color to be used for the message.
	 */
	@Override
	public void setMessage(String message, Color color) {
		//parent.setInfoText(message, color, true);
	}
}
