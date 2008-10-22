package fi.csc.microarray.databeans.features;

import java.util.List;

import fi.csc.microarray.MicroarrayException;

public class ModifiedFeature extends FeatureBase {

	private List<Feature> originals;

	public ModifiedFeature(List<Feature> originals) {
		this.originals = originals;
	}
	
	public Iterable<Float> asFloats() throws MicroarrayException {
		return null;
	}

	public Iterable<String> asStrings() throws MicroarrayException {
		return null;
	}

	public Table asTable() throws MicroarrayException {
		return null;
	}

	public boolean exists() {
		for (Feature original : originals) {
			if (!original.exists()) {
				return false;
			}
		}
		return true;
	}

	public String getName() {
		return null;
	}

}
