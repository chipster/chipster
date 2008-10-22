package fi.csc.microarray.databeans.features;

import java.util.LinkedList;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;

public class ConstantFloatFeature extends BasicFeature {

	private LinkedList<Float> floatList = new LinkedList<Float>();

	public ConstantFloatFeature(DataBean bean, FeatureProvider factory, float floatValue) {
		super(bean, factory);
		this.floatList.add(floatValue);
	}
	
	public Iterable<Float> asFloats() throws MicroarrayException {		
		return floatList;
	}

	public Iterable<String> asStrings() throws MicroarrayException {
		return null;
	}

	public Table asTable() throws MicroarrayException {
		return null;
	}
}
