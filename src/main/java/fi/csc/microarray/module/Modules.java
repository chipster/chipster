package fi.csc.microarray.module;

import java.util.LinkedList;

import fi.csc.microarray.databeans.DataManager;

public class Modules {
	
	
	private LinkedList<Module> modules = new LinkedList<Module>();

	public void addModule(Module module) {
		this.modules.add(module);
	}

	public void plugFeatures(DataManager manager) {
		for (Module module : this.modules) {
			module.plugFeatures(manager);
			module.plugModifiers(manager);
			module.plugContentTypes(manager);
		}
	}
}
