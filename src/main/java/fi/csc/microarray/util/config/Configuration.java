package fi.csc.microarray.util.config;

import java.util.HashMap;
import java.util.Map;

public class Configuration {
	
	private static final String PROPERTY_REFERENCE_POSTFIX = "}";
	private static final String PROPERTY_REFERENCE_PREFIX = "${";

	public Map<String, String[]> values = new HashMap<String, String[]>();
	public Map<String, Configuration> subModules = new HashMap<String, Configuration>();

	private boolean rewriteSystemPropertyReferences;
	
	public Configuration() {
		this(false);
	}
	
	public Configuration(boolean rewriteSystemPropertyReferences) {
		this.rewriteSystemPropertyReferences = rewriteSystemPropertyReferences;
	}

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
	
	public Configuration getModule(String name) {
		return subModules.get(name);
	}
	
	public void putValueInSubmodule(String moduleName, String name, String value) {
		Configuration module = getSubmoduleIfExists(moduleName);
		module.putValue(name, value);
	}

	public void putValue(String name, String value) {
		putValues(name, new String[] { value });
	}

	public void putValuesInSubmodule(String moduleName, String name, String[] values) {
		Configuration module = getSubmoduleIfExists(moduleName);
		module.putValues(name, values);
	}

	public void putValues(String name, String[] value) {
		if (rewriteSystemPropertyReferences) {
			for (int i = 0; i < value.length; i++) {
				if (value[i].startsWith(PROPERTY_REFERENCE_PREFIX)) {
					String[] split = value[i].split(PROPERTY_REFERENCE_POSTFIX);
					String systemProperty = System.getProperty(split[0].substring(2));
					value[i] = systemProperty + split[1];
				}
			}
		}
		values.put(name, value);
	}

	private Configuration getSubmoduleIfExists(String moduleName) {
		Configuration module = this;
		if (getModule(moduleName) != null) {
			module = getModule(moduleName);
		}
		return module;
	}
	

	public void createSubModule(String name) {
		subModules.put(name, new Configuration(this.rewriteSystemPropertyReferences));
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

}
