package fi.csc.microarray.client.operation.parameter;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.util.LinkedList;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.JComponent;
import javax.swing.JCheckBox;
import javax.swing.JPanel;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.parameter.EnumParameter.SelectionOption;

/**
 * A component that represents an EnumParameter that can have several
 * values selected at once. Each available options is shown as a
 * checkbox.
 * 
 * @author naktinis
 *
 */
@SuppressWarnings("serial")
public class MultipleSelectionInputComponent extends ParameterInputComponent
                                             implements ActionListener {
    
    // Logger for this class
    private static final Logger logger = Logger.
        getLogger(SingleSelectionInputComponent.class);
    
    EnumParameter param;
    List<JCheckBox> checkboxes = new LinkedList<JCheckBox>();
    JPanel checkboxPanel = new JPanel();

    protected MultipleSelectionInputComponent(EnumParameter param, 
                                              ParameterPanel parameterPanel) {
        super(parameterPanel);
        this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
        this.param = param;
        
        // Generate a list of checkboxes representing available choices
        boolean selected;
        for (Object choice : param.getOptions()) {
            SelectionOption option = (SelectionOption) choice;
            
            // Select default values
            selected = false; 
            if (param.getSelectedOptions().contains(option)) {
                selected = true;
            }
            
            // Create the component
            JCheckBox checkbox = new JCheckBox(option.toString(), selected);
            checkbox.setName(option.getValue());
            checkbox.addActionListener(this);
            checkbox.addFocusListener(this);
            
            // Add checkbox to lists
            checkboxes.add(checkbox);
            this.add(checkbox);
        }
        this.getComponent(0);
            
    }
    
    /**
     * @return a list of names of selected checkboxes. 
     */
    public List<String> getSelectedNames() {
        List<String> selectionList = new LinkedList<String>();
        for (JCheckBox checkbox : checkboxes) {
            if (checkbox.isSelected()) {
                selectionList.add(checkbox.getName());
            }
        }
        return selectionList;
    }

    @Override
    public Parameter getParameter() {
        return param;
    }

    @Override
    public JComponent getParameterComponent() {
        return this;
    }

    @Override
    public boolean inputIsValid() {
        /**
         * Currently we assume that any of the given options is
         * automatically valid as an input. We only check the
         * number of chosen values.
         */
        int count = getSelectedNames().size();
        return (count >= param.getMaxCount()) && (count <= param.getMinCount());
    }

    public void focusGained(FocusEvent arg0) {
        getParentPanel().setMessage(param.getDescription(), Color.black);
    }

    public void actionPerformed(ActionEvent arg0) {
        param.setValue(getSelectedNames());
        logger.debug("set selection to " + getSelectedNames());
        getParentPanel().setMessage(param.getDescription(), Color.black);
    }

}
