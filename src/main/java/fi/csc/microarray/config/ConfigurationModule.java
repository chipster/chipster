package fi.csc.microarray.config;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

public class ConfigurationModule {
	
	private static final String PROPERTY_REFERENCE_POSTFIX = "}";
	private static final String PROPERTY_REFERENCE_PREFIX = "${";

	public static final String[] VALUE_MUST_BE_SET = new String[] {};
	
	public Map<String, String[]> values = new HashMap<String, String[]>();
	public Map<String, ConfigurationModule> subModules = new HashMap<String, ConfigurationModule>();

	public String[] getValues(String name) {
		return values.get(name);
	}

	public String getValue(String name) {
		if (values.get(name) == null) {
			return null;
			
		} else if (values.get(name).length > 1) {
			throw new IllegalArgumentException("multiple values found for " + name);
			
		} else {
			return values.get(name)[0];
		}		
	}
	
	public ConfigurationModule getModule(String name) {
		return subModules.get(name);
	}
	
	public void putValueInSubmodule(String moduleName, String name, String value) {
		ConfigurationModule module = getModule(moduleName);
		module.putValue(name, value);
	}

	public void putValue(String name, String value) {
		putValues(name, new String[] { value });
	}

	public void putValuesInSubmodule(String moduleName, String name, String[] values) {
		ConfigurationModule module = getModule(moduleName);
		module.putValues(name, values);
	}

	public void putValues(String name, String[] value) {
		for (int i = 0; i < value.length; i++) {
			if (value[i].startsWith(PROPERTY_REFERENCE_PREFIX)) {
				String[] split = value[i].split(PROPERTY_REFERENCE_POSTFIX);
				String systemProperty = System.getProperty(split[0].substring(2));
				value[i] = systemProperty + split[1];
			}
		}
		values.put(name, value);
	}

	public void createSubModule(String name) {
		subModules.put(name, new ConfigurationModule());
	}

	public void addValue(String name, String value) {
		String[] oldValues = values.get(name);
		
		if (oldValues == null) {
			putValue(name, value);
			
		} else {
			String[] newValues = new String[oldValues.length +1];
			for (int i = 0; i < oldValues.length; i++) {
				newValues[i] = oldValues[i];
			}
			newValues[newValues.length-1] = value;
			putValues(name, newValues);
		}
	}

	public boolean hasSubModule(String name) {
		return subModules.containsKey(name);
	}

	public List<String> findMissingValues() {
		LinkedList<String> missingValues = new LinkedList<String>();
		
		// check this module
		for (String key : values.keySet()) {
			// check if some values are not yet set
			if (values.get(key) == VALUE_MUST_BE_SET) {
				missingValues.add(key);
			}
		}
		
		// recurse into other modules
		for (ConfigurationModule module : subModules.values()) {
			missingValues.addAll(module.findMissingValues());
		}
		
		return missingValues;
	}

	public void removeValues(String name) {
		values.remove(name);		
	}

}
