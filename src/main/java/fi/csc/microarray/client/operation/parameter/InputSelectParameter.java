package fi.csc.microarray.client.operation.parameter;

import java.util.List;

import fi.csc.microarray.client.operation.Operation.DataBinding;

public class InputSelectParameter extends DataSelectionParameter {

	protected InputSelectParameter(String name, String description, String initValue) {
		super(name, description, initValue);
	}

	public void setDataBindings(List<DataBinding> bindings) {
		SelectionOption[] inputs = new SelectionOption[bindings.size()];
		int initIndex = 0;
		for (int i = 0; i < bindings.size(); i++) {
			inputs[i] = new SelectionOption(bindings.get(i).getData().getName(), bindings.get(i).getName());
			if (initValue != null && initValue.equals(bindings.get(i).getName())) {
				initIndex = i;
			}
		}
		setOptions(inputs, initIndex);	
	}
}
