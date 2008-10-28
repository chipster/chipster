package fi.csc.microarray;

import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.Modules;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.module.stats.StatsModule;

public class ModulesForTesting {
	
	private static Module[] modulesForTesting = new Module[] {
			new BasicModule(), 
			new StatsModule(), 
			new MicroarrayModule()
	};
	
	public static Modules getModulesForTesting() {
		Modules modules = new Modules();
		for (Module module : modulesForTesting) {
			modules.addModule(module);
		}
		return modules;
	}

}
