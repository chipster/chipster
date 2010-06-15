package fi.csc.microarray.databeans.features;

import fi.csc.microarray.databeans.Dataset;

public class BoolTrueFeature extends ConstantStringFeature  {

	public BoolTrueFeature(Dataset bean, FeatureProvider factory) {
		super(bean, factory, "true");
	}
	
}
