package fi.csc.microarray.module;

import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.module.stats.StatsModule;

public class DefaultModules {

	public static Modules getDefaultModules() {
		Modules modules = new Modules();
		modules.addModule(new BasicModule());
		modules.addModule(new StatsModule());
		modules.addModule(new MicroarrayModule());
		return modules;
	}	

}
