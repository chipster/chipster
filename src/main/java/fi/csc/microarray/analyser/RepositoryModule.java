package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Category;
import fi.csc.microarray.util.XmlUtil;

/**
 * <p>One module in repository, corresponds to one (modulename)-module.xml file.</p>
 * 
 * <p>Access is strictly synchronised, because all operations
 * may lead to module or script updates if files on the disk have changed. To avoid
 * deadlocking, dependencies must be kept one way: RepositoryModule never calls ToolRepository and
 * ToolDescription never calls RepositoryModule.</p>
 *
 * @author Taavi Hupponen, Aleksi Kallio 
 */
public class RepositoryModule {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(RepositoryModule.class);
	
	private LinkedList<Category> categories = new LinkedList<Category>();
	private LinkedHashSet<ToolDescription> descriptions = new LinkedHashSet<ToolDescription>();
	private LinkedHashSet<String> supportedDescriptions = new LinkedHashSet<String>();
	private LinkedHashSet<String> visibleDescriptions = new LinkedHashSet<String>();

	private File moduleDir;
	private File moduleFile;
	private long moduleFileTimestamp;
	private HashMap<String, ToolRuntime> runtimes;

	private String summary = null;
	private String moduleName = null;
	
	public RepositoryModule(File moduleDir, File moduleFile, HashMap<String, ToolRuntime> runtimes) throws ParserConfigurationException, FileNotFoundException, SAXException, IOException {
		this.moduleFile = moduleFile;
		this.moduleDir = moduleDir;
		this.runtimes = runtimes;
		load();
	}
	
	public synchronized ModuleDescriptionMessage getModuleDescriptionMessage() {
		
		// Check if up-to-date
		reloadModuleIfNeeded();		
		
		// Construct description message using the current state 
		ModuleDescriptionMessage msg = new ModuleDescriptionMessage(moduleName);
		for (Category category : categories) {
			msg.addCategory(category);
		}
		
		return msg;
	}

	public synchronized ToolDescription getDescription(String id) {

		// Check if up-to-date
		reloadModuleIfNeeded();		
		
		// Find description
		ToolDescription desc = findDescription(id);
		
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
		return findDescription(id); 
	}

	private ToolDescription findDescription(String id) {
		
		// Always iterate over descriptions, because they can change their ID on the fly
		for (ToolDescription description : descriptions) {
			if (id.equals(description.getID())) {
				return description;
			}
		}
		
		// Matching description was not found
		return null;
	}

	public synchronized boolean isSupportedDescription(String id) {

		// Check if up-to-date
		reloadModuleIfNeeded();		

		return supportedDescriptions.contains(id);
	}

	public synchronized File getModuleDir() {
		return moduleDir;
	}

	public synchronized String getSummary() throws FileNotFoundException, SAXException, IOException, ParserConfigurationException {
		
		// Check if up-to-date
		reloadModuleIfNeeded();		

		return summary;
	}

	private synchronized void updateDescription(ToolDescription desc) {
	    // FIXME params should not be empty
	    HashMap<String, String> params = new HashMap<String, String>();
		ToolDescription newDescription;
		try {
			newDescription = desc.getHandler().handle(moduleDir, desc.getToolFile().getName(), params);
		} catch (AnalysisException e) {
			// update failed, continue using the old one
			return;
		}
		if (newDescription != null) {
			newDescription.setUpdatedSinceStartup();
			
			// name (id) of the tool has not changed
			if (desc.getID().equals(newDescription.getID())) {
				
				// replace the old description with the same name
				descriptions.add(newDescription);
				if (supportedDescriptions.contains(desc.getID())) {
					supportedDescriptions.add(newDescription.getID());
				}
				if (visibleDescriptions.contains(desc.getID())) {
					visibleDescriptions.add(newDescription.getID());
				}
			} 

			// name (id) of the tool has changed
			else {
				logger.warn("ID of the tool has changed, registering it with old and new name");
				if (findDescription(newDescription.getID()) != null){
					logger.warn("descriptions already contain a tool with the new name, ignoring the new tool");
					return;
				} 
				// add the tool with the new name
				descriptions.add(newDescription);
				if (supportedDescriptions.contains(desc.getID())) {
					supportedDescriptions.add(newDescription.getID());
				}
				if (visibleDescriptions.contains(desc.getID())) {
					visibleDescriptions.add(newDescription.getID());
				}
			}
		}
	}

	/**
	 * Parses a module file and loads all tools listed in it.
	 * 
	 * Disabled tools are not available in this analyser instance,
	 * but are available in some other instances running, so the
	 * client still has to be informed about them.
	 * 
	 * @param moduleDir directory where module specification file and runtime specific tool directories are located in
	 * @param moduleFile module specification file
	 * @return summary of tool counts
	 * 
	 * @throws FileNotFoundException
	 * @throws SAXException
	 * @throws IOException
	 * @throws ParserConfigurationException
	 */
	private synchronized void load() throws FileNotFoundException, SAXException, IOException, ParserConfigurationException {

		
		// Update timestamp, so that even if loading fails, we do not try to load again before the file has been changed (hopefully fixed)
		moduleFileTimestamp = moduleFile.lastModified();
		
		// Load module description file
		Document document = XmlUtil.parseReader(new FileReader(moduleFile));
		Element moduleElement = (Element)document.getElementsByTagName("module").item(0);
		
		// Load and check module name 
		this.moduleName  = moduleElement.getAttribute("name");
		if (moduleName.isEmpty()) {
			this.summary  = "not loading a module without a name";
			logger.warn(summary);
			return;
		}
		
		// Initialise stats
	    int totalCount = 0;
		int successfullyLoadedCount = 0;
	    int hiddenCount = 0;
	    int disabledCount = 0;
		
		// Load categories and tools in them
		for (Element categoryElement: XmlUtil.getChildElements(moduleElement, "category")) {

			// Category name
			String categoryName = categoryElement.getAttribute("name");
			if (categoryName.isEmpty()) {
				logger.warn("not loading a category without a name");
				continue;
			}
			
			// Enabled or disabled status
			boolean categoryDisabled = categoryElement.getAttribute("disabled").equals("true");
			if (categoryDisabled) {
				logger.info("not loading category " + categoryName + ": disabled");
				continue;
			}
			
			// GUI color
			String categoryColor = categoryElement.getAttribute("color");
			if (categoryColor.isEmpty()) {
				logger.warn("not loading category " + categoryName + ": no color");
				continue;
			}
			
		    // Category visibility
		    boolean categoryHidden = Boolean.valueOf(categoryElement.getAttribute("hidden"));

		    // Create and register the category
		    Category category = new Category(categoryName, categoryColor, categoryHidden);
		    categories.add(category);
		    
		    // Load tools and add them to category
		    for (Element toolElement: XmlUtil.getChildElements(categoryElement, "tool")) {
		    	totalCount++;

		    	// Resource (script file name, Java class name, ...)
		    	Element resourceElement = XmlUtil.getChildElement(toolElement, "resource");
		    	if (resourceElement == null) {
		    		logger.warn("not loading a tool without resource element");
		    		continue;
		    	}
		    	String resource = resourceElement.getTextContent().trim();
		    	if (resource == null || resource.isEmpty()) {
		    		logger.warn("not loading a tool with empty resource element");
		    		continue;
		    	}

		    	// Tool visibility 
		    	boolean toolDisabled = toolElement.getAttribute("disabled").equals("true");
		    	boolean toolHidden = categoryHidden; // currently tools can be hidden only if their category is hidden

		    	// Tool runtime
		    	String runtimeName = toolElement.getAttribute("runtime");
		    	ToolRuntime runtime = runtimes.get(runtimeName); // FIXME depends on repository!!!
		    	if (runtime == null) {
		    		logger.warn("not loading " + resource + ": runtime " + runtimeName + " not found");
		    		continue;
		    	}

		    	// Tool parameters
		    	boolean parametersOk = true;
		    	HashMap<String, String> parameters = new HashMap<String, String>();
		    	for (Element parameterElement : XmlUtil.getChildElements(toolElement, "parameter")) {
		    		String parameterName = XmlUtil.getChildElement(parameterElement, "name").getTextContent().trim();
		    		if (parameterName == null || parameterName.isEmpty()) {
		    			logger.warn("parameter without a name");
		    			parametersOk = false;
		    			break;
		    		}

		    		String parameterValue = XmlUtil.getChildElement(parameterElement, "value").getTextContent().trim();
		    		if (parameterValue == null) {
		    			logger.warn("parameter without a value");
		    			parametersOk = false;
		    			break;
		    		}

		    		// This parameter is ok
		    		parameters.put(parameterName, parameterValue);
		    	}
		    	if (!parametersOk) {
		    		logger.warn("not loading " + resource + ": parameters not ok");
		    		continue;
		    	}
		    	
		    	// Create the analysis description
		    	ToolDescription description;
		    	try {
		                description = runtime.getHandler().handle(moduleDir, resource, parameters);
		                
		                // check for duplicates in this module
				        ToolDescription previousDescription = getDescription(description.getID());
			    	    if (previousDescription != null) {
			    	        logger.warn("not loading " + resource + ": tool with the same ID already exists in this module");
			    	        continue;
			    	    }
		    	    
		    	} catch (AnalysisException e) {
		    		logger.warn("loading " + resource + " failed, could not create description", e);
		    		continue;
		    	}
		    	
		    	// Register the tool
		    	descriptions.add(description);
		    	successfullyLoadedCount++;

		    	// Set disabled if needed
		    	String disabledStatus = "";
		    	if (!runtime.isDisabled() && !toolDisabled) {
		    		// Not disabled, add to supported descriptions list
		    		supportedDescriptions.add(description.getID());
		    		
		    	} else {
		    		disabledStatus = " DISABLED";
		    		disabledCount++;
		    	}

	    		// Add to category, which gets sent to the client
	    		category.addTool(description.getSADL(), description.getHelpURL());
	    		
                // Set hidden if needed    		
                String hiddenStatus = "";
		    	if (toolHidden) {
		    		hiddenStatus = " HIDDEN";
		    		hiddenCount++;
		    	}

		    	logger.info("loaded " + description.getID() + " " + description.getDisplayName() + " " +
		    			description.getToolFile() + disabledStatus + hiddenStatus);
		    }

		}

		// Update summary
		this.summary = "loaded " + moduleName + " " + successfullyLoadedCount + "/" + totalCount +
		" tools, " + disabledCount + " disabled, " + hiddenCount + " hidden";
		logger.info(summary);
	}
	
	private synchronized void reloadModuleIfNeeded() {
		if (moduleFile.lastModified() > moduleFileTimestamp) {
			
			// Clean everything
			this.categories.clear();
			this.descriptions.clear();
			this.moduleName = null;
			this.summary = null;
			this.supportedDescriptions.clear();
			this.visibleDescriptions.clear();
			
			// Reload
			try {
				load();
				
			} catch (Exception e) {
				logger.error("reloading module " + moduleName + " failed", e);
			}
		}
	}

}