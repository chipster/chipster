package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.util.LinkedList;
import java.util.List;

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;
import javax.swing.UIManager;

import org.jdesktop.swingx.JXTaskPane;
import org.jdesktop.swingx.JXTaskPaneContainer;
import org.jdesktop.swingx.VerticalLayout;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationPanel;
import fi.csc.microarray.client.operation.OperationDefinition.InputDefinition;
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
@SuppressWarnings("serial")
public class ToolParameterPanel extends ParameterPanel {
	private JScrollPane scroller;

	private static final int LEFT_MARGIN = 10;
	private static final int TOP_MARGIN = 2;

	private OperationPanel parent;
	
	/**
	 * Creates a new ParameterPanel for the given operation.
	 * 
	 * @param operation The operation which is to be started from this panel.
	 * @param client The client that is to be informed when a job is done.
	 * @throws MicroarrayException 
	 */
	public ToolParameterPanel(Operation operation, OperationPanel parent) throws MicroarrayException {

	    super(operation, new BorderLayout());
		this.parent = parent;
		
		// Configure style of parameter panel
        UIManager.put("TaskPaneContainer.background",
            UIManager.getColor("TaskPane.background"));
        UIManager.put("TaskPane.borderColor",
            UIManager.getColor("TaskPane.background"));
		
		// Create a collapsible pane container
        JXTaskPaneContainer paneContainer = new JXTaskPaneContainer();
        paneContainer.setBorder(null);
        JXTaskPane pane;
        JPanel paramPane;
        GridBagConstraints con;

        // Remove vertical gap
        VerticalLayout verticalLayout = new VerticalLayout();
        verticalLayout.setGap(0);
        paneContainer.setLayout(verticalLayout);
        paneContainer.setBackground(this.getBackground());
	    
        
        // Divide parameters into required and optional
        List<Parameter> requiredParameters = new LinkedList<Parameter>();
        List<Parameter> optionalParameters = new LinkedList<Parameter>();
        for (Parameter param : operation.getParameters()) {
            if (param.isOptional()) {
                optionalParameters.add(param);
            } else {
                requiredParameters.add(param);
            }
        }

        // Parameters
        paramPane = new JPanel(new GridBagLayout());
        con = prepareBagConstraints(); 
        
        // Required parameters
        if (requiredParameters.size() > 0) {
    		for (Parameter param : requiredParameters) {
    			ParameterInputComponent component = createInputComponent(param);
                JLabel label = component.getLabel();
                label.setFont(label.getFont().deriveFont(label.getFont().getStyle() ^ Font.BOLD));
    			addParameter(paramPane, component, label, con);
    		}
    		
            // Add required parameters to the collapsible pane
	        paneContainer.add(paramPane);
		}
        
        // Optional parameters
        if (optionalParameters.size() > 0) {
            for (Parameter param : optionalParameters) {
                ParameterInputComponent component = createInputComponent(param);
                addParameter(paramPane, component, component.getLabel(), con);
            }
            
            // Add optional parameters to the collapsible pane
            paneContainer.add(paramPane);
        }
        
        // Input file mappings
        pane = new JXTaskPane();
        pane.setTitle("Input datasets");
        pane.setCollapsed(false);

        pane.getContentPane().setBackground(this.getBackground());
        pane.setBackground(this.getBackground());
        
        // Grid layout for component/label pairs
        paramPane = new JPanel(new GridBagLayout());
        paramPane.setBackground(this.getBackground());
        paramPane.setBorder(BorderFactory.createEmptyBorder());

        con = prepareBagConstraints();
        
        List<InputFileComponent> inputComponents = new LinkedList<InputFileComponent>();
        
        // Only show input mappings in parameter panel when necessary
        if (operation.getDefinition().getInputs().size() > 1) {
            // Operation has some inputs
            for (InputDefinition input : operation.getDefinition().getInputs()) {
                InputFileComponent inputComponent = new InputFileComponent(input, operation);
                inputComponent.setListener(inputComponent.new InputFileComponentListener(inputComponents));
                inputComponents.add(inputComponent);
                
                addParameter(paramPane, inputComponent, inputComponent.getLabel(), con);
            }
                  
            // Add the inputs to the collapsable pannel
            pane.add(paramPane);
            paneContainer.add(pane);
            
        }

		scroller = new JScrollPane(paneContainer);
		scroller.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		
		
		scroller.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
	    
		this.add(scroller, BorderLayout.CENTER);
	}

	/**
	 * Routine for adding component/label pair to a panel.
	 * 
	 * @param panel Panel that will contain the component.
	 * @param component Control that will be added.
	 * @param label JLabel object defining.
	 * @param con Constraint object that defines Control's position.
	 */
	private void addParameter(JPanel panel, JComponent component, JLabel label,
	                          GridBagConstraints con) {       
        con.gridx = 0;
        con.gridy++;
        con.insets.top = TOP_MARGIN;
        con.insets.left = LEFT_MARGIN;
        con.fill = GridBagConstraints.HORIZONTAL;
        panel.add(label, con);
        con.gridx = 1;
        con.anchor = GridBagConstraints.EAST;
        con.fill = GridBagConstraints.NONE;
        panel.add(component, con);
	}
	
	/**
	 * Utility routine for preparing GridBagConstraints.
	 * 
	 * @param con
	 * @return initialized GridBagConstraint.
	 */
	private GridBagConstraints prepareBagConstraints() {
	    GridBagConstraints con = new GridBagConstraints();
        con.gridx = 0; con.gridy = 0;
        con.gridwidth = 1;
        con.weightx = 1.0; con.weighty = 0;
        con.anchor = GridBagConstraints.WEST;
        return con;
	}
	
	
	/**
	 * Sets the message of the (bottom left) infobox of this parameter panel.
	 * 
	 * @param message The text to be set.
	 * @param color The font color to be used for the message.
	 */
	@Override
	public void setMessage(String message, Color color) {
		parent.setInfoText(message, color, true);
	}
}
