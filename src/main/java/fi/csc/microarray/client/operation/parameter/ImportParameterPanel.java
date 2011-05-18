package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.JPanel;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.ToolPanel;
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
 * @author Janne Käki, Aleksi Kallio, Petri Klemelä
 *
 */
@SuppressWarnings("serial")
public class ImportParameterPanel extends ParameterPanel {

	private static final int LEFT_MARGIN = 10;
	private static final int TOP_MARGIN = 2;
	
	
	
	/**
	 * Creates a new ParameterPanel for the given operation.
	 * 
	 * @param operation The operation which is to be started from this panel.
	 * @param client The client that is to be informed when a job is done.
	 * @throws MicroarrayException 
	 */
	public ImportParameterPanel(Operation operation, ToolPanel parent) throws MicroarrayException {
		super(operation, new BorderLayout());
	
		JPanel paramPane = new JPanel(new GridBagLayout());
		GridBagConstraints con = new GridBagConstraints();
		
		con.gridx = 0; con.gridy = 0;
		con.gridwidth = 1;
		con.weightx = 1.0; con.weighty = 0;
		con.anchor = GridBagConstraints.WEST;
		
		
		for (Parameter param : operation.getParameters()) {
			ParameterInputComponent component = createInputComponent(param);
					
			con.gridx = 0;
			con.gridy++;
			con.insets.top = TOP_MARGIN;
			con.insets.left = LEFT_MARGIN;
			con.fill = GridBagConstraints.NONE;
			paramPane.add(component.getLabel(), con);
			con.gridx = 1;
			con.anchor = GridBagConstraints.WEST;
			con.fill = GridBagConstraints.NONE;
			paramPane.add(component, con);
		}
		
		con.weighty = 1;
		con.weightx = 0;
		con.gridx = 0;		
		con.gridy++;
		con.gridwidth = 2;
		con.fill = GridBagConstraints.BOTH;
		paramPane.add(new JPanel(),con);
		
		this.add(paramPane, BorderLayout.WEST);

	
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
