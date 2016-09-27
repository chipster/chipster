package fi.csc.microarray.client.operation.parameter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
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
				InputParameters inputParameters, List<DataBean> options, DataBean selected) {
			super(input);
			
	        // Prepare the combo box
	        choiceBox = new SteppedComboBox(options.toArray());
			Dimension preferredSize = choiceBox.getPreferredSize();
			choiceBox.setMinimumSize(ParameterInputComponent.PREFERRED_SIZE);	
			choiceBox.setPreferredSize(ParameterInputComponent.PREFERRED_SIZE);
			choiceBox.setPopupWidth(preferredSize.width);
			choiceBox.setBackground(Color.white);
	        
	        // Set selected bean for this combo box, if null -> no selection
			choiceBox.setSelectedItem(selected);
	        
	        this.add(choiceBox, BorderLayout.CENTER);
	        
	        this.initListener(inputParameters);
		}
		
	    private void initListener(InputParameters inputParameters) {
	        choiceBox.addItemListener(new ItemListener() {
				@Override
				public void itemStateChanged(ItemEvent e) {
					if (e.getStateChange() == ItemEvent.DESELECTED) {
		                /* deselection happens only when the data is selected in some 
		                 * other comboBox and in this case that comboBox will set the
		                 * new data bindings
		                 */
		            } else {
		            	inputParameters.bind((DataBean)e.getItem(), e.getSource());
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
		private List<DataBean> options;

		public MultiInput(InputDefinition input, Operation operation,
				InputParameters inputParameters, List<DataBean> options, List<DataBean> selected) {
					
			super(input);
			
			this.options = options;
			
	        list = new JList<DataBean>(options.toArray(new DataBean[0]));
			//Dimension preferredSize = list.getPreferredSize();
	        list.setMinimumSize(new Dimension((int) ParameterInputComponent.PREFERRED_SIZE.getWidth(), list.getMinimumSize().height));
	        list.setBorder(new LineBorder(Color.gray));
			//choiceBox.setPopupWidth(preferredSize.width);
	        //choiceBox.setBackground(Color.white);
	        
	        // Set selected bean for this combo box, if null -> no selection
	        list.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
	        	      
	        for (int i = 0; i < selected.size(); i++) {
	        	if (selected.contains(options.get(i))) {
	        		list.addSelectionInterval(i, i);
	        	}
	        }
	       
	        this.add(list, BorderLayout.CENTER);
	        
    		this.initListener(inputParameters);
		}
		
		private void initListener(InputParameters inputParameters) {
			list.addListSelectionListener(new ListSelectionListener() {
				@Override
				public void valueChanged(ListSelectionEvent e) {
					
					for (DataBean data : list.getSelectedValuesList()) {
						inputParameters.bind(data, e.getSource());
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
    		int index = options.indexOf(selected);
			list.removeSelectionInterval(index, index);
		}
	}
	
    
    @SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(SingleSelectionInputComponent.class);
	private InputDefinition input;
    
    public InputFileComponent(InputDefinition input) {
        super(new BorderLayout());
        
        // Set name and operation object
        this.input = input;
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
}