package fi.csc.microarray.databeans;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

import javax.swing.Icon;

import org.apache.log4j.Logger;

import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProvider;
import fi.csc.microarray.databeans.features.Modifier;
import fi.csc.microarray.exception.MicroarrayException;

public abstract class DataManagerBase implements DataManager {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(DataManagerBase.class);

	private Map<String, FeatureProvider> factories = new HashMap<String, FeatureProvider>();
	private Map<String, Modifier> modifiers = new HashMap<String, Modifier>();
	
	/** MIME types for the DataBeans */
	private Map<String, ContentType> contentTypes = new HashMap<String, ContentType>();
	
	/** Mapping file extensions to content types */
	private Map<String, String> extensionMap = new HashMap<String, String>();
	
	private LinkedList<DataChangeListener> listeners = new LinkedList<DataChangeListener>();
	
	private boolean eventsEnabled = false;

	public void addDataChangeListener(DataChangeListener listener) {
		logger.debug("adding DataChangeListener: " + listener);
		if (listener == null) {
			throw new IllegalArgumentException("listener cannot be null");
		}
		listeners.add(listener);
	}

	public void setEventsEnabled(boolean enabled) {
		this.eventsEnabled = enabled;		
	}

	public void close() {
		// do nothing
	}
	
	public void dispatchEventIfVisible(DataChangeEvent event) {
		if (event.getDataItem().getParent() != null) {
			dispatchEvent(event);
		}
	}
	
	public void dispatchEvent(DataChangeEvent event) {
		if (eventsEnabled) {
			// dispatch events only for connected datas
			for (DataChangeListener listener : listeners) {
				if (listener == null) {
					logger.error("One of the DataChangeListeners listeners was null.");
				} else {
					logger.debug("Notifying DataChangeListener " + listener.toString());
				}
				listener.dataChanged(event);
			}
		}
	}
	
	public DataBean createDataBean(String name, DataFolder parent, DataBean source) throws MicroarrayException {
		return createDataBean(name, parent, wrapSource(source));
	}
	
	public static DataBean[] wrapSource(DataBean source) {
		DataBean[] sources = null;
		
		if (source != null) {
			sources = new DataBean[1];
			sources[0] = source;
		} else {
			sources = new DataBean[0];
		}
		
		return sources;
		
	}
	
	public void plugContentType(String mimeType, boolean supported, boolean binary, String description, Icon icon, String... extensions) {
		// create the content type
		contentTypes.put(mimeType, new ContentType(mimeType, supported, binary, description, icon, extensions));
		
		
		// add extensions to search map
		for (String extension: extensions) {
			extensionMap.put(extension, mimeType);
		}
	}
	
	public void plugModifier(String name, Modifier modifier) {
		modifiers.put(name, modifier);
	}
	
	public void plugFeatureFactory(String name, FeatureProvider plugin) {
		logger.debug("plugged " + plugin.getClass().getSimpleName() + " at " + name);
		plugin.setName(name);
		factories.put(name, plugin);
	}
	
	
	public Modifier fetchModifier(String modifierName) {	
		return modifiers.get(modifierName);
	}
	
	public Feature fetchFeature(String featureName, DataBean bean) {
		String bestMatch = null; 		
		for (String feature : factories.keySet()) {
			if (featureName.startsWith(feature)) {
				if (bestMatch == null || feature.length() > bestMatch.length()) {
					// current best match
					bestMatch = feature;
				}
			}
		}
		FeatureProvider factory = factories.get(bestMatch);
		if (factory == null) {
			throw new RuntimeException("no feature factory plugged in for \"" + featureName + "\" (total of " + factories.size() + " factories plugged)");
		}
		logger.debug("best match for " + featureName + " was " + (factory != null ? factory.getName() : factory));
		String namePostfix = getNamePostfix(featureName, factory.getName());		
		return factory.createFeature(namePostfix, bean);
	}

	/**
	 * Guess the MIME content type using the filename.
	 * 
	 * For now, simply use the extension to figure out the mime type.
	 * 
	 * Types are plugged at ApplicationConstants.
	 * 
	 */
	public ContentType guessContentType(String name) {
		
		ContentType type = null;
		if(name.lastIndexOf(".") != -1){
			String extension = name.substring(name.lastIndexOf(".") + 1, name.length()).toLowerCase();
			String typeName = extensionMap.get(extension);
			if (typeName != null) {
				type = contentTypes.get(typeName);
			}
		} 
		
		if (type == null) {
			type = contentTypes.get("application/octet-stream");
		}
		return type;
	}
	
	
	public ContentType guessContentType(File file) {
		return guessContentType(file.getName());
	}
	
	public ContentType getContentType(String typeName) {
		return contentTypes.get(typeName);
	}
	
	public DataItem findDataItem(String name) {
		return findDataItem(name, getRootFolder());
	}
	
	private DataItem findDataItem(String name, DataItem root) {
		DataItem matchingItem = null;
		
		// root item matches
		if (root.getName().equals(name)) {
			return root;
		} 
		
		// root is a folder, search children
		else if (root instanceof DataFolder) {
			for (DataItem child: ((DataFolder)root).getChildren()) {
				matchingItem = findDataItem(name, child);
				if (matchingItem != null) {
					return matchingItem;
				}
			}
		} 
		
		// no match found
		return null;
	}
	
	
	
	private String getNamePostfix(String featureName, String factoryName) {
		if (factoryName.length() > featureName.length()) {
			return "";
		} else {
			String npf = featureName.substring(factoryName.length());
			if (npf.startsWith("/")) {
				return npf.substring(1);
			} else {
				return npf;
			}
		}
	}
	
	
}
