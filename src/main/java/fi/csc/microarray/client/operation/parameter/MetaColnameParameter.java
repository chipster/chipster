package fi.csc.microarray.client.operation.parameter;

import java.util.List;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;

public class MetaColnameParameter extends ColnameParameter {

	protected MetaColnameParameter(String name, String description, String initValue) {
		super(name, description, initValue);
	}

	public void setDataBindings(List<DataBinding> bindings) throws MicroarrayException {
		DataBean metaData = null;
		for (DataBinding binding : bindings) {
			if (binding.getInputType() == ChipsterInputTypes.PHENODATA) {
				metaData = binding.getData();
				break;
			}
		}
		loadOptionsFromColumnNames(metaData);
	}
}
