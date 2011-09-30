package fi.csc.microarray.module;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.module.basic.BasicModule;

public class ModuleManager {
	
	private LinkedList<Module> modules = new LinkedList<Module>();
	private Module primaryModule = null;
	
	public ModuleManager() throws InstantiationException, IllegalAccessException, ClassNotFoundException {
		this(null);
	}
	
	public ModuleManager(String primaryModuleClass) throws InstantiationException, IllegalAccessException, ClassNotFoundException {

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
	
	public void plugAll(DataManager manager, Session session) {
		for (Module module : this.modules) {
			module.plugFeatures(manager);
			module.plugModifiers(manager);
			module.plugContentTypes(manager);
			module.plugTypeTags(manager);
			if (session != null) {
				session.getVisualisations().addVisualisationMethods(module.getVisualisationMethods());
			}
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
