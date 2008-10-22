package fi.csc.microarray.databeans.features.stat;

import java.util.List;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.features.CalculatingIterable;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.ModifiedFeature;
import fi.csc.microarray.databeans.features.CalculatingIterable.CalcOperation;

public class CalcModifier {

	protected static class CalcModifierFeature extends ModifiedFeature {
	
			private Feature original;
			private CalcOperation operation;
			
			protected CalcModifierFeature(List<Feature> inputs, CalcOperation operation) {
				super(inputs);
				if (inputs.size() != 1) {
					throw new IllegalArgumentException("calc modifier must have 1 parameter");
				}
				this.original = inputs.get(0);
				this.operation = operation;
			}
	
			public Iterable<Float> asFloats() throws MicroarrayException {			
				return new CalculatingIterable(original.asFloats(), null, operation);
			}
		}

	private CalcModifierFeature output;
	private final CalcOperation operation;

	public CalcModifier(CalcOperation operation) {
		this.operation = operation;
	}

	public Feature getOutput() {
		return output;
	}

	public void setInputs(List<Feature> inputs) {
		this.output = new CalcModifierFeature(inputs, operation);
	
	}
}