package fi.csc.microarray.module;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.module.basic.BasicModule;

public class Modules {
	
	private LinkedList<Module> modules = new LinkedList<Module>();
	private Module primaryModule = null;
	
	public Modules() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		this(null);
	}
	
	public Modules(String primaryModuleClass) throws InstantiationException, IllegalAccessException, ClassNotFoundException {

		// load default module
		modules.add(new BasicModule());

		// load primary module, if given
		if (primaryModuleClass != null) {
			Module primaryModule = (Module)Class.forName(primaryModuleClass).newInstance();
			this.primaryModule = primaryModule;
			modules.add(primaryModule);
			
		} else {
			this.primaryModule = modules.getLast(); // use the last default module as primary
		}
	}
	
	public void plugFeatures(DataManager manager) {
		for (Module module : this.modules) {
			module.plugFeatures(manager);
			module.plugModifiers(manager);
			module.plugContentTypes(manager);
			module.plugTypeTags(manager);
		}
	}

	public void addModule(Module module) {
		this.modules.add(module);
	}
	
	public List<Module> getModules() {
		return modules;
	}

	public Module getPrimaryModule() {
		return primaryModule;
	}
}
