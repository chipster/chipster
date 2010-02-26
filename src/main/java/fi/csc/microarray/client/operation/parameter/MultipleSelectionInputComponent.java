package fi.csc.microarray.client.operation.parameter;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JCheckBox;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.parameter.EnumParameter.SelectionOption;

/*
 * TODO: Create a JPanel which encapsulates the checkboxes. This
 *       JPanel should be returned by getParameterComponent method.
 */

@SuppressWarnings("serial")
public class MultipleSelectionInputComponent extends ParameterInputComponent
                                             implements ActionListener {
    
    // Logger for this class
    private static final Logger logger = Logger.
        getLogger(SingleSelectionInputComponent.class);
    
    EnumParameter param;
    List<JCheckBox> checkboxes = new LinkedList<JCheckBox>();

    protected MultipleSelectionInputComponent(EnumParameter param, 
                                              ParameterPanel parameterPanel) {
        super(parameterPanel);
        this.param = param;
        
        // Generate a list of checkboxes representing available choices
        for (Object choice : param.getOptions()) {
            SelectionOption option = (SelectionOption) choice;
            checkboxes.add(new JCheckBox(option.toString()));
        }
            
    }

    @Override
    public Parameter getParameter() {
        return param;
    }

    @Override
    public JComponent getParameterComponent() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public boolean inputIsValid() {
        // TODO Auto-generated method stub
        return false;
    }

    public void focusGained(FocusEvent arg0) {
        // TODO Auto-generated method stub
        
    }

    public void actionPerformed(ActionEvent arg0) {
        // TODO Auto-generated method stub
        
    }

}
