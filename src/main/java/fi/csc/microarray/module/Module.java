package fi.csc.microarray.module;

import fi.csc.microarray.databeans.DataManager;

public interface Module {

	public void plugFeatures(DataManager manager);
	public void plugModifiers(DataManager manager);
	public void plugContentTypes(DataManager manager);
}
