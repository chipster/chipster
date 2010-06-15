package fi.csc.microarray.client.operation.parameter;

import java.util.List;

import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;

public class ColnameParameter extends DataSelectionParameter {

	protected ColnameParameter(String id, String displayName, String description, String initValue) {
		super(id, displayName, description, initValue);
	}

	public void setDataBindings(List<DataBinding> bindings) throws MicroarrayException {
		Dataset data = null;
		for (DataBinding binding : bindings) {
			if (binding.getInputType() != ChipsterInputTypes.PHENODATA) {
				data = binding.getData();
				break;
			}
		}
		loadOptionsFromColumnNames(data);
	}
}
