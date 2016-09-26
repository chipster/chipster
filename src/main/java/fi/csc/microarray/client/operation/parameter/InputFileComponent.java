package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.ListSelectionModel;
import javax.swing.border.LineBorder;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

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
public abstract class InputFileComponent extends JPanel {
	
	
	public static class SingleInput extends InputFileComponent{
		
		private SteppedComboBox choiceBox;

		public SingleInput(InputDefinition input, Operation operation,
				List<InputFileComponent> components) {
			super(input, operation, components);
			
	        // Prepare the combo box
	        choiceBox = new SteppedComboBox(dataBeans.toArray());
	        choiceBox.setEnabled(enabled);
			Dimension preferredSize = choiceBox.getPreferredSize();
			choiceBox.setMinimumSize(ParameterInputComponent.PREFERRED_SIZE);	
			choiceBox.setPreferredSize(ParameterInputComponent.PREFERRED_SIZE);
			choiceBox.setPopupWidth(preferredSize.width);
			choiceBox.setBackground(Color.white);
	        
	        // Set selected bean for this combo box, if null -> no selection
	        if (!currentBeans.isEmpty()) {
	        	choiceBox.setSelectedItem(currentBeans.get(0));
	        } else {
	        	choiceBox.setSelectedItem(null);
	        }
	        
	        this.add(choiceBox, BorderLayout.CENTER);
	        
	        this.initListener(this);
		}
		
	    private void initListener(InputFileComponent parent) {
	        choiceBox.addItemListener(new ItemListener() {
				@Override
				public void itemStateChanged(ItemEvent e) {
					if (e.getStateChange() == ItemEvent.DESELECTED) {
		                /* deselection happens only when the data is selected in some 
		                 * other comboBox and in this case that comboBox will set the
		                 * new data bindings
		                 */
		            } else {
		            	parent.bind((DataBean)e.getItem(), e.getSource());
		            }
				}
			});
		}

		public JComponent getComponent() {
	        return choiceBox;
	    }
	    
		public List<DataBean> getSelectedItems() {
	    	List<DataBean> items = new ArrayList<>();
	    	for (Object dataBean : choiceBox.getSelectedObjects()) {
	    		items.add((DataBean)dataBean);
	    	}
			return items;
	    }
	    

		public void removeSelected(Object selected) {
			if (selected.equals(choiceBox.getSelectedItem())) {
				choiceBox.setSelectedItem(null);
			}
		}
	}
	
	public static class MultiInput extends InputFileComponent {
		
		private JList<DataBean> list;

		public MultiInput(InputDefinition input, Operation operation,
				List<InputFileComponent> components) {
			super(input, operation, components);
			
	        list = new JList<DataBean>(dataBeans.toArray(new DataBean[0]));
	        list.setEnabled(enabled);
			//Dimension preferredSize = list.getPreferredSize();
	        list.setMinimumSize(new Dimension((int) ParameterInputComponent.PREFERRED_SIZE.getWidth(), list.getMinimumSize().height));
	        list.setBorder(new LineBorder(Color.gray));
			//choiceBox.setPopupWidth(preferredSize.width);
	        //choiceBox.setBackground(Color.white);
	        
	        // Set selected bean for this combo box, if null -> no selection
	        list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
	        	      
	        for (int i = 0; i < dataBeans.size(); i++) {
	        	if (currentBeans.contains(dataBeans.get(i))) {
	        		list.addSelectionInterval(i, i);
	        	}
	        }
	       
	        this.add(list, BorderLayout.CENTER);
	        
    		this.initListener(this);
		}
		
		private void initListener(InputFileComponent parent) {
			list.addListSelectionListener(new ListSelectionListener() {
				@Override
				public void valueChanged(ListSelectionEvent e) {
					for (int i = e.getFirstIndex(); i <= e.getLastIndex(); i++) {
						// skip deselection events
						if (list.isSelectedIndex(i)) {
							parent.bind(dataBeans.get(i), e.getSource());
						}
					}
				}
			});
		}
		
	    public JComponent getComponent() {
	        return list;
	    }
	    
	    public List<DataBean> getSelectedItems() {
    		return list.getSelectedValuesList();	
	    }
	    

		public void removeSelected(Object selected) {
    		int index = dataBeans.indexOf(selected);
			list.removeSelectionInterval(index, index);
		}
	}
	
    
    private static final Logger logger = Logger
        .getLogger(SingleSelectionInputComponent.class);
    
    private InputDefinition input;
    private Operation operation;
    
	protected List<DataBean> dataBeans;

	private List<InputFileComponent> components;

	protected ArrayList<DataBean> currentBeans;

	protected boolean enabled;

    public InputFileComponent(InputDefinition input, Operation operation, List<InputFileComponent> components) {
        super(new BorderLayout());
        
        // Set name and operation object
        this.input = input;
        this.operation = operation;
        this.components = components;
        
        // Check current bindings and generate choices
        List<DataBinding> bindings = operation.getBindings();
        HashMap<String, String> bindingMap = new HashMap<String, String>();
        DataBean[] dataBeans = new DataBean[0];
        currentBeans = new ArrayList<DataBean>();
        enabled = false;
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
                    currentBeans.add(binding.getData());
                }
            }
        }
        
        this.dataBeans = Arrays.asList(dataBeans);
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
    abstract public JComponent getComponent();
    abstract public List<DataBean> getSelectedItems();
	abstract public void removeSelected(Object selected);
        
    /**
     * @return parameter name for this component.
     */
    public JLabel getLabel() {
        return new JLabel(input.getDisplayName());
    }
    
    private void bind(DataBean selected, Object sourceComponent) {

        logger.debug("Selected input dataset: " + selected );

        operation.clearBindings();
        LinkedList<DataBinding> newBindings = new LinkedList<DataBinding>();
        for (InputFileComponent component : components) {
            // Make sure no other input control has the same value
            if (component.getSelectedItems().contains(selected) 
            		&& sourceComponent != component.getComponent()) {
                component.removeSelected(selected);;
            }
            
            // Rebind input datasets                    
            List<DataBean> selectedBeans = component.getSelectedItems();
            
            InputDefinition input = component.getInput();
            
            input.resetMulti();
            
            for (DataBean bean : selectedBeans) {
            		newBindings.add(new DataBinding(bean, input.getID(), input.getType()));
            		input.nextMulti();
            }
        }
        operation.setBindings(newBindings);
    }
}