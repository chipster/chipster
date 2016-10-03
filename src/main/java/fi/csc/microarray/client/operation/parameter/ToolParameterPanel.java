package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;

import javax.swing.BorderFactory;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ScrollPaneConstants;

import org.jdesktop.swingx.JXPanel;
import org.jdesktop.swingx.JXTaskPane;
import org.jdesktop.swingx.JXTaskPaneContainer;
import org.jdesktop.swingx.VerticalLayout;
import org.jdesktop.swingx.painter.MattePainter;

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
public class ToolParameterPanel extends ParameterPanel {
	private JScrollPane scroller;

	private static final int LEFT_MARGIN = 10;
	private static final int TOP_MARGIN = 2;

	private ToolPanel parent;
	
	/**
	 * Creates a new ParameterPanel for the given operation.
	 * 
	 * @param operation The operation which is to be started from this panel.
	 * @param client The client that is to be informed when a job is done.
	 * @throws MicroarrayException
	 */ 
	public ToolParameterPanel(Operation operation, ToolPanel parent) throws MicroarrayException {

	    super(operation, new BorderLayout());
		this.parent = parent;
		
                
		// Create a collapsible pane container
        JXTaskPaneContainer paneContainer = new JXTaskPaneContainer();
        JXTaskPane pane;
        JPanel paramPane;
        GridBagConstraints con;

        // Remove vertical gap
        VerticalLayout verticalLayout = new VerticalLayout();
        verticalLayout.setGap(0);
        paneContainer.setLayout(verticalLayout);
	    
        // Parameters
        paramPane = new JPanel(new GridBagLayout());
        con = prepareBagConstraints(); 
        for (Parameter param : operation.getParameters()) {
            if (param.isOptional()) {
                ParameterInputComponent component = createInputComponent(param);
                addParameter(paramPane, component, component.getLabel(), con);
            } else {
    			ParameterInputComponent component = createInputComponent(param);
                JLabel label = component.getLabel();
                label.setFont(label.getFont().deriveFont(label.getFont().getStyle() ^ Font.BOLD));
    			addParameter(paramPane, component, label, con);

            }
        }
        
        paneContainer.add(paramPane);
        
        // Input file mappings
        pane = new JXTaskPane();
        pane.setTitle("Input datasets");
        pane.setCollapsed(false);

        
        // Grid layout for component/label pairs
        paramPane = new JPanel(new GridBagLayout());

        con = prepareBagConstraints();        
        
        // Only show input mappings in parameter panel when necessary
        if (operation.getDefinition().getInputs().size() > 1) {
            // Operation has some inputs
        	
        	InputParameters inputs = new InputParameters(operation);
        	
            for (InputFileComponent inputComponent : inputs.getComponents()) {
                addParameter(paramPane, inputComponent, inputComponent.getLabel(), con);
            }
                  
            // Add the inputs to the collabsiple panel
            pane.add(paramPane);
            paneContainer.add(pane);
            
        }

		scroller = new JScrollPane(paneContainer);
		scroller.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		
		// without this background painter fills all extra space with its default background (blue in windows)
		paneContainer.setBackgroundPainter(new MattePainter(this.getBackground()));
		pane.getContentPane().setBackground(this.getBackground());
		// change JXPanel default border which draws some small quirks around the panel
		((JXPanel)pane.getContentPane()).setBorder(BorderFactory.createEmptyBorder(5, 10, 10, 10));
		
		paneContainer.setBorder(BorderFactory.createEmptyBorder());

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
		parent.updateSuitability();
	}
}
