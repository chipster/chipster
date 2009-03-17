package fi.csc.microarray.config;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;


public class ConfigurationModule {
	
	public Map<String, ConfigurationEntry> values = new HashMap<String, ConfigurationEntry>();
	public Map<String, ConfigurationModule> subModules = new HashMap<String, ConfigurationModule>();

	public ConfigurationEntry getEntry(String name) {
		return values.get(name);
	}
	
	public void addEntry(ConfigurationEntry entry) {
		values.put(entry.getName(), entry);
	}
	
	public ConfigurationModule getModule(String name) {
		return subModules.get(name);
	}

	public void createSubModule(String name) {
		subModules.put(name, new ConfigurationModule());
	}

	public boolean hasSubModule(String name) {
		return subModules.containsKey(name);
	}

	public List<String> findMissingValues() {
		LinkedList<String> missingValues = new LinkedList<String>();
		
		// check this module
		for (String key : values.keySet()) {
			// check if some values are not yet set
			if (values.get(key).mustBeSet()) {
				missingValues.add(key);
			}
		}
		
		// recurse into other modules
		for (ConfigurationModule module : subModules.values()) {
			missingValues.addAll(module.findMissingValues());
		}
		
		return missingValues;
	}
}
