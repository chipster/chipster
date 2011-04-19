package fi.csc.microarray.analyser;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Category;

public class RepositoryModule {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(RepositoryModule.class);
	
	private LinkedHashMap<String, AnalysisDescription> descriptions = new LinkedHashMap<String, AnalysisDescription>();
	private LinkedHashSet<String> supportedDescriptions = new LinkedHashSet<String>();
	private LinkedHashSet<String> visibleDescriptions = new LinkedHashSet<String>();

	private ModuleDescriptionMessage moduleDescriptionMessage;

	private File moduleFile;
	private File moduleDir;
	
	public RepositoryModule(File moduleDir, File moduleFile, String moduleName) throws ParserConfigurationException {
		this.moduleFile = moduleFile;
		this.moduleDir = moduleDir;
		this.moduleDescriptionMessage = new ModuleDescriptionMessage(moduleName);
	}
	
	public void updateDescription(AnalysisDescription desc) throws AnalysisException {
	    // FIXME params should not be empty
	    HashMap<String, String> params = new HashMap<String, String>();
		AnalysisDescription newDescription = desc.getHandler().handle(null, desc.getToolFile().toString(), params);
		if (newDescription != null) {
			newDescription.setUpdatedSinceStartup();
			
			// name (id) of the tool has not changed
			if (desc.getID().equals(newDescription.getID())) {
				
				// replace the old description with the same name
				descriptions.put(newDescription.getID(), newDescription);
				if (supportedDescriptions.contains(desc.getID())) {
					supportedDescriptions.add(newDescription.getID());
				}
				if (visibleDescriptions.contains(desc.getID())) {
					visibleDescriptions.add(newDescription.getID());
				}
			} 

			// name (id) of the tool has changed
			else {
				logger.warn("name of the tool has changed after loading from custom-scripts, keeping both old and new");
				if (descriptions.containsKey(newDescription.getID())){
					logger.warn("descriptions already contains a tool with the new name, ignoring custom-scripts");
					return;
				} 
				// add the tool with the new name
				descriptions.put(newDescription.getID(), newDescription);
				if (supportedDescriptions.contains(desc.getID())) {
					supportedDescriptions.add(newDescription.getID());
				}
				if (visibleDescriptions.contains(desc.getID())) {
					visibleDescriptions.add(newDescription.getID());
				}
			}
		}
	}

	public void addDescriptionCategory(Category category) {
		getModuleDescriptionMessage().addCategory(category);		
	}

	public ModuleDescriptionMessage getModuleDescriptionMessage() {
		return moduleDescriptionMessage;
	}

	public AnalysisDescription getDescription(String id) {
		return descriptions.get(id);
	}

	public void putDescription(String id, AnalysisDescription description) {
		descriptions.put(id, description);
	}

	public boolean hasDescription(String id) {
		return descriptions.containsKey(id);
	}

	public void markSupportedDescription(String id) {
		supportedDescriptions.add(id);
		
	}

	public boolean isSupportedDescription(String id) {
		return supportedDescriptions.contains(id);
	}

	
}