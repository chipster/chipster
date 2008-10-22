package fi.csc.microarray.databeans.features;

import java.util.LinkedList;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;

public class ConstantStringFeature extends BasicFeature {

	private LinkedList<String> stringList = new LinkedList<String>();

	public ConstantStringFeature(DataBean bean, FeatureProvider factory, String string) {
		super(bean, factory);
		this.stringList.add(string);
	}
	
	public Iterable<Float> asFloats() throws MicroarrayException {
		return null;
	}

	public Iterable<String> asStrings() throws MicroarrayException {
		return stringList;
	}

	public Table asTable() throws MicroarrayException {
		return null;
	}
}
