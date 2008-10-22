package fi.csc.microarray.databeans.features;

import fi.csc.microarray.databeans.DataBean;

public class BoolTrueFeature extends ConstantStringFeature  {

	public BoolTrueFeature(DataBean bean, FeatureProvider factory) {
		super(bean, factory, "true");
	}
	
}
