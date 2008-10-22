package fi.csc.microarray.wizard;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class WizardParameterBundle {

	private HashMap<String, Map<String, List<String>>> mappedStrings = new HashMap<String, Map<String,List<String>>>();
	private HashMap<String, String> strings = new HashMap<String, String>();
	
	public void add(String name, Map<String, List<String>> strings) {
		if (mappedStrings.containsKey(name)) {
			throw new IllegalStateException(name + " already added");
		}
		mappedStrings.put(name, strings);
	}

	public void add(String name, String string) {
		if (strings.containsKey(name)) {
			throw new IllegalStateException(name + " already added");
		}
		strings.put(name, string);
	}

	public Map<String, List<String>> getMappedStrings(String name) {
		return mappedStrings.get(name);
	}

	public String getString(String name) {
		return strings.get(name);
	}

}
