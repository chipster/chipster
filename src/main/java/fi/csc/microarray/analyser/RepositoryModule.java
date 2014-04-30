
package fi.csc.microarray.analyser;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.analyser.ToolDescription.InputDescription;
import fi.csc.microarray.analyser.ToolDescription.OutputDescription;
import fi.csc.microarray.config.DirectoryLayout;
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
	
	private LinkedList<CategoryInModule> categories = new LinkedList<CategoryInModule>();
	private LinkedHashMap<String, ToolDescription> descriptions = new LinkedHashMap<String, ToolDescription>();
	private LinkedHashSet<String> supportedDescriptions = new LinkedHashSet<String>();
	private LinkedHashSet<String> visibleDescriptions = new LinkedHashSet<String>();

	private File moduleDir;
	private File moduleFile;
	private long moduleFileTimestamp;
	private HashMap<String, ToolRuntime> runtimes;

	private String summary = null;
	private String moduleName = null;
	
    public static class CategoryInModule {
        private String name;
        private String color;
        private Boolean hidden;
        private List<ToolDescription> tools = new LinkedList<ToolDescription>();
        
        public CategoryInModule(String name, String color, Boolean hidden) {
            this.name = name;
            this.color = color;
            this.hidden = hidden;
        }
        
        public String getName() {
            return name;
        }
        
        public String getColor() {
            return color;
        }
        
        public Boolean isHidden() {
            return hidden;
        }
        
        public void addTool(ToolDescription tool) {
            tools.add(tool);
        }
        
        public List<ToolDescription> getTools() {
            return tools;
        }
    }

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
		
		for (CategoryInModule categoryInModule : categories) {
			Category category = new Category(categoryInModule.getName(), categoryInModule.getColor(), categoryInModule.isHidden());
			
			for (ToolDescription tool : categoryInModule.getTools()) {
				ToolDescription freshTool = reloadToolIfNeeded(tool); // check that we are using up-to-date description
				category.addTool(freshTool.getSADL(), freshTool.getHelpURL());
			}
			msg.addCategory(category);
		}
		
		return msg;
	}

	public synchronized ToolDescription getDescription(String id) {

		// Check if up-to-date
		reloadModuleIfNeeded();		
		
		// Find description
		ToolDescription desc = descriptions.get(id);
		
		// Return null if nothing is found
		if (desc == null) {
			return null;
		}

		// Check if description needs to be updated
		reloadToolIfNeeded(desc);
		
		// Return the possibly updated description
		return descriptions.get(id); 
	}

	private ToolDescription reloadToolIfNeeded(ToolDescription oldDesc) {
		
		if (oldDesc != null && !oldDesc.isUptodate()) {
			ToolDescription newDesc = updateDescription(oldDesc);
			logger.info("updated tool: " + oldDesc.getID());
			return newDesc;
			
		} else {
			return oldDesc;	
		}
		
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

	private synchronized ToolDescription updateDescription(ToolDescription oldDescription) {
	    // FIXME params should not be empty
	    HashMap<String, String> params = new HashMap<String, String>();
		ToolDescription newDescription;
		try {
			newDescription = oldDescription.getHandler().handle(oldDescription.getModuleDir(), oldDescription.getToolFile().getName(), params);
			
		} catch (AnalysisException e) {
			// update failed, continue using the old one
			return oldDescription;
		}
		
		if (newDescription == null) {
			return oldDescription;
		}

		newDescription.setUpdatedSinceStartup();

		// name (id) of the tool has not changed
		if (oldDescription.getID().equals(newDescription.getID())) {

			// replace the old description with the same name
			descriptions.put(newDescription.getID(), newDescription);

			// FIXME what's the point?
			if (supportedDescriptions.contains(oldDescription.getID())) {
				supportedDescriptions.add(newDescription.getID());
			}
			if (visibleDescriptions.contains(oldDescription.getID())) {
				visibleDescriptions.add(newDescription.getID());
			}
		} 

		// name (id) of the tool has changed
		// FIXME maybe safer not to proceed
		else {
			logger.warn("name of the tool was changed, keeping both old and new");
			if (descriptions.containsKey(newDescription.getID())){
				logger.warn("descriptions already contains a tool with the new name, ignoring reload");
				return oldDescription;
			} 
			// add the tool with the new name
			descriptions.put(newDescription.getID(), newDescription);
			if (supportedDescriptions.contains(oldDescription.getID())) {
				supportedDescriptions.add(newDescription.getID());
			}
			if (visibleDescriptions.contains(oldDescription.getID())) {
				visibleDescriptions.add(newDescription.getID());
			}
		}

		return newDescription;
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
		    CategoryInModule category = new CategoryInModule(categoryName, categoryColor, categoryHidden);
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
		    	ToolRuntime runtime = runtimes.get(runtimeName);
		    	if (runtime == null) {
		    		logger.warn("not loading " + resource + ": runtime " + runtimeName + " not found");
		    		continue;
		    	}

		    	// Tool module
		    	File toolModuleDir;
		    	String nonDefaultModuleName = toolElement.getAttribute("module");
		    	if (nonDefaultModuleName != null && !nonDefaultModuleName.equals("")) {
		    		toolModuleDir = new File(DirectoryLayout.getInstance().getModulesDir(), nonDefaultModuleName);
		    	} else {
		    		toolModuleDir = moduleDir;
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
		            description = runtime.getHandler().handle(toolModuleDir, resource, parameters);
		    	    
		    	} catch (AnalysisException e) {
		    		logger.warn("loading " + resource + " failed, could not create description", e);
		    		continue;
		    	}
		    	
		    	// Check that filenames are unique. Overwriting input files is a bad idea when the input file is
		    	// only a symlink to the original file 		
		    	boolean filenamesOk = true;
		    	HashSet<String> uniqueNames = new HashSet<>();
		    	ArrayList<String> allNames = new ArrayList<>();
		    	
		    	for (InputDescription input : description.getInputFiles()) {
		    		allNames.add(input.getFileName());
		    	}
		    	
		    	for (OutputDescription output : description.getOutputFiles()) {
		    		allNames.add(output.getFileName().getID());
		    	}
		    	
		    	for (String name : allNames) {
		    		if (name == null) {
		    			// name is null for file sets
		    			continue;
		    		}
		    		if (uniqueNames.contains(name)) {
		    			logger.warn("filename " + name + " isn't unique");
		    			filenamesOk = false;
		    		} else {
		    			uniqueNames.add(name);
		    		}
		    	}
		    	
		    	if (!filenamesOk) {
		    		logger.warn("not loading " + resource + ": non-unique filename(s)");
		    		continue;
		    	}
		    	
		    	// Register the tool, override existing
		    	descriptions.put(description.getID(), description);
		    	successfullyLoadedCount++;

		    	// Set disabled if needed, override existing
		    	String disabledStatus = "";
		    	if (!runtime.isDisabled() && !toolDisabled) {
		    		// Not disabled, add to supported descriptions list
		    		supportedDescriptions.add(description.getID());
		    		
		    	} else {
		   			supportedDescriptions.remove(description.getID());	
		    		disabledStatus = " DISABLED";
	    			disabledCount++;
		    	}

	    		// Add to category, which gets sent to the client
	    		category.addTool(description);
	    		
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