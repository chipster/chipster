package fi.csc.microarray.client.operation.parameter;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationDefinition.InputDefinition;
import fi.csc.microarray.databeans.DataBean;

public class InputParameters {
	
	private List<InputFileComponent> inputComponents = new LinkedList<InputFileComponent>();
	private Operation operation;
	
	public InputParameters(Operation operation) {
		this.operation = operation;

		// create component for each input
		for (InputDefinition input : operation.getDefinition().getInputs()) {
			
			List<DataBean> options = new ArrayList<>();
			List<DataBean> selected = new ArrayList<>();
			
			// collect selected and available databeans
			for (DataBinding binding : operation.getBindings()) {
				if (OperationDefinition.doBackwardsCompatibleTypeCheck(input.getType(), binding.getData())) {
					options.add(binding.getData());
				}
				
				// Find the selected item
				if (input.idMatches(binding.getName())) {
					selected.add(binding.getData());
				}
			}
			
			InputFileComponent inputComponent;
			if (input.isMulti()) {
				inputComponent = new InputFileComponent.MultiInput(input, operation, this, options, selected);
			} else {
				// optional parameters may not have binding
				DataBean singleSelected = selected.isEmpty() ? null : selected.get(0);
				inputComponent = new InputFileComponent.SingleInput(input, operation, this, options, singleSelected);
			}
			inputComponents.add(inputComponent);
		}
	}
	
	public List<InputFileComponent> getComponents() {
		return inputComponents;
	}	
	
    public void bind(DataBean selected, Object sourceComponent) {

        // Make sure no other input control has the same value
        for (InputFileComponent component : inputComponents) {
            if (component.getSelectedItems().contains(selected) 
            		&& sourceComponent != component.getComponent()) {
                component.removeSelected(selected);;
            }
        }
        
        LinkedList<DataBinding> newBindings = new LinkedList<DataBinding>();
        
        // collect bindings from the components                    
        for (InputFileComponent component : inputComponents) {
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
