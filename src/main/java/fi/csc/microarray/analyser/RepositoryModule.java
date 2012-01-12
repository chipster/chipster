package fi.csc.microarray.analyser;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Category;

/**
 * <p>One module in repository, corresponds to one (modulename)-module.xml file.</p>
 * 
 * <p>Access is strictly synchronised, because all operations
 * may lead to module or script updates if files on the disk have changed. To avoid
 * deadlocking, dependencies must be kept one way: RepositoryModule never calls ToolRepository and
 * AnalysisDescription never calls RepositoryModule.</p>
 *
 * @author Taavi Hupponen, Aleksi Kallio 
 */
public class RepositoryModule {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(RepositoryModule.class);
	
	private LinkedHashMap<String, AnalysisDescription> descriptions = new LinkedHashMap<String, AnalysisDescription>();
	private LinkedHashSet<String> supportedDescriptions = new LinkedHashSet<String>();
	private LinkedHashSet<String> visibleDescriptions = new LinkedHashSet<String>();

	private ModuleDescriptionMessage moduleDescriptionMessage;

	private File moduleDir;
	
	public RepositoryModule(File moduleDir, File moduleFile, String moduleName) throws ParserConfigurationException {
		// moduleFile is not used yet, but needed for update checks that are implemented later
		this.moduleDir = moduleDir;
		this.moduleDescriptionMessage = new ModuleDescriptionMessage(moduleName);
	}
	
	public void updateDescription(AnalysisDescription desc) throws AnalysisException {
	    // FIXME params should not be empty
	    HashMap<String, String> params = new HashMap<String, String>();
		AnalysisDescription newDescription = desc.getHandler().handle(this, desc.getToolFile().getName(), params);
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
		return moduleDescriptionMessage; // FIXME remove this and construct on fly
	}

	public AnalysisDescription getDescription(String id) throws AnalysisException {
		
		// Find description
		AnalysisDescription desc = descriptions.get(id);
		
		// Return null if nothing is found
		if (desc == null) {
			return null;
		}

		// Check if description needs to be updated
		if (desc != null && !desc.isUptodate()) {
			updateDescription(desc);
			logger.info("updated tool: " + desc.getID());
		}
		
		// Return the possibly updated description
		return descriptions.get(id); 
	}

	public void putDescription(String id, AnalysisDescription description) {
		descriptions.put(id, description);
	}

	public void markSupportedDescription(String id) {
		supportedDescriptions.add(id);
		
	}

	public boolean isSupportedDescription(String id) {
		return supportedDescriptions.contains(id);
	}

	public void setModuleDir(File moduleDir) {
		this.moduleDir = moduleDir;
	}

	public File getModuleDir() {
		return moduleDir;
	}

	
}