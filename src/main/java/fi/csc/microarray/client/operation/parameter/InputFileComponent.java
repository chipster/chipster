package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.OperationDefinition.InputDefinition;
import fi.csc.microarray.databeans.DataBean;

/**
 * User can remap input files to appropriate input
 * parameters.
 * 
 * @author naktinis
 *
 */
@SuppressWarnings("serial")
public class InputFileComponent extends JPanel {
    
    private static final Logger logger = Logger
        .getLogger(SingleSelectionInputComponent.class);
    
    private InputDefinition input;
    private SteppedComboBox choiceBox;
    private Operation operation;

    public InputFileComponent(InputDefinition input, Operation operation) {
        super(new BorderLayout());
        
        // Set name and operation object
        this.input = input;
        this.operation = operation;
        
        // Check current bindings and generate choices
        List<DataBinding> bindings = operation.getBindings();
        HashMap<String, String> bindingMap = new HashMap<String, String>();
        DataBean[] dataBeans = new DataBean[0];
        DataBean currentBean = null;
        Boolean enabled = false;
        Integer index = 0;
        if (!bindings.isEmpty()) {
            
            // User has already selected the input files
            enabled = true;
            dataBeans = new DataBean[bindings.size()];
            for (DataBinding binding : bindings) {
                dataBeans[index++] = binding.getData();
                bindingMap.put(binding.getName(), binding.getData().getName());
                
                // Find the selected item
                if (input.idMatches(binding.getName())) {
                    currentBean = binding.getData();
                }
            }
        }
        
        // Prepare the combo box
        choiceBox = new SteppedComboBox(dataBeans);
        choiceBox.setEnabled(enabled);
		Dimension preferredSize = choiceBox.getPreferredSize();
        choiceBox.setPreferredSize(ParameterInputComponent.PREFERRED_SIZE);
		choiceBox.setPopupWidth(preferredSize.width);
        choiceBox.setBackground(Color.white);
        
        // Set selected bean for this combo box, if null -> no selection
        choiceBox.setSelectedItem(currentBean);
        
        this.add(choiceBox, BorderLayout.CENTER);
    }
    
    /**
     * @return input definition object.
     */
    public InputDefinition getInput() {
        return input;
    }
    
    /**
     * @return checkbox component.
     */
    public JComboBox getChoiceBox() {
        return choiceBox;
    }
    
    /**
     * Sets action listener for this component.
     */
    public void setListener(InputFileComponentListener listener) {
        choiceBox.addItemListener(listener);
    }
    
    /**
     * @return parameter name for this component.
     */
    public JLabel getLabel() {
        return new JLabel(input.getDisplayName());
    }
    
    /**
     * Listener that is responsible for changes made to input
     * file selection components.
     */
    public class InputFileComponentListener implements ItemListener {
        
        private Object selected;
        private List<InputFileComponent> components;
        
        /**
         * Set a list that contains InputFileComponent objects
         * for this operation. Each InputFileComponent needs to
         * manipulate the full list. Should include itself also.
         */
        public InputFileComponentListener(List<InputFileComponent> components) {
            this.components = components;
        }

        /**
         * Basically does two things:
         * <ol>
         * <li> Makes sure that no other InputFileComponent has the same value.
         * <li> Rebinds inputs in the Operation object.
         * </ol>
         */
        public void itemStateChanged(ItemEvent e) {
            if (e.getStateChange() == ItemEvent.DESELECTED) {
                /* deselection happens only when the data is selected in some 
                 * other comboBox and in this case that comboBox will set the
                 * new data bindings
                 */
            } else {
                selected = e.getItem();

                logger.debug("Selected input dataset: " + selected );

                operation.clearBindings();
                LinkedList<DataBinding> newBindings = new LinkedList<DataBinding>();
                for (InputFileComponent component : components) {
                    // Make sure no other input control has the same value
                    if (component.getChoiceBox().getSelectedItem() != null 
                    		&& component.getChoiceBox().getSelectedItem().equals(selected) 
                    		&& e.getSource() != component.getChoiceBox()) {
                        component.getChoiceBox().setSelectedItem(null);
                    }
                    
                    // Rebind input datasets                    
                    DataBean selectedBean = (DataBean) component.getChoiceBox().getSelectedItem();
                    if (selectedBean != null) {
                    	newBindings.add(new DataBinding(selectedBean,
                    			component.getInput().getID(),
                    			component.getInput().getType()));
                    }
                }
                operation.setBindings(newBindings);
            }
        }
    }
}