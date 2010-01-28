package fi.csc.microarray.databeans.features;

import java.util.Iterator;

import fi.csc.microarray.exception.MicroarrayException;

public abstract class FeatureBase implements Feature {

	public String asString() throws MicroarrayException {
		Iterator<String> strings = asStrings().iterator();
		String string = strings.next();
		if (strings.hasNext()) {
			throw new IllegalStateException("more than one string was returned");
		}
		return string;
	}
	
	public Float asFloat() throws MicroarrayException {
		Iterator<Float> floats = asFloats().iterator();
		Float floatValue = floats.next();
		if (floats.hasNext()) {
			throw new IllegalStateException("more than one float was returned");
		}
		return floatValue;
	}

}