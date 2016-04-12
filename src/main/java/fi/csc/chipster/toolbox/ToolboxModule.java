
package fi.csc.chipster.toolbox;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import de.schlichtherle.truezip.file.TFileReader;
import fi.csc.chipster.toolbox.SADLTool.ParsedScript;
import fi.csc.chipster.toolbox.toolpartsparser.ToolPartsParser;
import fi.csc.chipster.util.StringUtils;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Category;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;
import fi.csc.microarray.util.XmlUtil;

/**
 * One module in toolbox , corresponds to one (modulename)-module.xml file.
 * 
 */
public class ToolboxModule {

	private static final Logger logger = Logger.getLogger(ToolboxModule.class);
	
	private LinkedList<ToolboxCategory> categories = new LinkedList<ToolboxCategory>();
	private LinkedHashMap<String, ToolboxTool> tools = new LinkedHashMap<String, ToolboxTool>();

	private File moduleDir;
	private File moduleFile;

	private String summary = null;
	private String moduleName = null;
	
    public static class ToolboxCategory {
        private String name;
        private String color;
        private Boolean hidden;
        private List<ToolboxTool> tools = new LinkedList<ToolboxTool>();
        
        public ToolboxCategory(String name, String color, Boolean hidden) {
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
        
        public void addTool(ToolboxTool tool) {
            tools.add(tool);
        }
        
        public List<ToolboxTool> getTools() {
        	return tools;
        }
        
    }

	public ToolboxModule(File moduleDir, File moduleFile) throws ParserConfigurationException, FileNotFoundException, SAXException, IOException {
		this.moduleFile = moduleFile;
		this.moduleDir = moduleDir;
		load();
	}
	
	public ModuleDescriptionMessage getModuleDescriptionMessage() {
		
		// Construct description message using the current state 
		ModuleDescriptionMessage msg = new ModuleDescriptionMessage(moduleName);
		
		for (ToolboxCategory toolboxCategory : categories) {
			Category category = new Category(toolboxCategory.getName(), toolboxCategory.getColor(), toolboxCategory.isHidden());
			
			for (ToolboxTool tool : toolboxCategory.getTools()) {
				
				// help url not supported (or used) at the moment
				category.addTool(tool.getSadlString(), null);
			}
			msg.addCategory(category);
		}
		
		return msg;
	}

	public ToolboxTool getTool(String id) {
		return tools.get(id);
	}

	public Collection<ToolboxTool> getAll() {
		return tools.values();
	}
	

	public String getSummary() {
		return summary;
	}

	
	/**
	 * Parses a module file and loads all tools listed in it.
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
	private void load() throws FileNotFoundException, SAXException, IOException, ParserConfigurationException {
		
		// Load module description file
		Document document = XmlUtil.parseReader(new TFileReader(moduleFile));
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
		    ToolboxCategory category = new ToolboxCategory(categoryName, categoryColor, categoryHidden);
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
		    	boolean toolHidden = categoryHidden; // currently tools can be hidden only if their category is hidden

		    	// Tool runtime
		    	String runtimeName = toolElement.getAttribute("runtime");
		    	if (runtimeName == null || runtimeName.isEmpty()) {
		    		logger.warn("not loading " + resource + ": runtime " + runtimeName + " is null or empty");
		    		continue;
		    	}

		    	// Tool module
		    	File toolModuleDir;
		    	String nonDefaultModuleName = toolElement.getAttribute("module");
		    	if (nonDefaultModuleName != null && !nonDefaultModuleName.equals("")) {
		    		toolModuleDir = new File(moduleDir.getParentFile(), nonDefaultModuleName);
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

		    	// get parser
	    		ToolPartsParser partsParser = Toolbox.getToolPartsParser(runtimeName);
	    		if (partsParser == null) {
	    			logger.warn("not loading " + resource + ": no parser for runtime " + runtimeName);
	    			continue;
	    		}

		    	// parse script parts (sadl / source)
	    		ParsedScript parsedScript;
	    		try {
		    		parsedScript = partsParser.parse(toolModuleDir, resource);
		    	} catch (Exception e) {
		    		logger.warn("loading " + resource + " failed, parsing parts failed", e);
		    		continue;
		    	}

		    	// check SADL by parsing it, also get tool id		
	    		SADLDescription sadlDescription;
	    		try {
	    			sadlDescription= new ChipsterSADLParser().parse(parsedScript.SADL);
		    	} catch (ParseException e) {
		    		logger.warn("loading " + resource + " failed, parsing sadl failed", e);
		    		continue;
		    	}
		    	
		    	String toolId = sadlDescription.getName().getID();
	    		
		    	// Check that filenames are unique. Overwriting input files is a bad idea when the input file is
		    	// only a symlink to the original file 		
		    	boolean filenamesOk = true;
		    	HashSet<String> uniqueNames = new HashSet<>();
		    	ArrayList<String> allNames = new ArrayList<>();
		    	
		    	for (Input input : sadlDescription.getInputs()) {
		    		allNames.add(input.getName().getID());
		    	}
		    	
		    	for (Output output : sadlDescription.getOutputs()) {
		    		allNames.add(output.getName().getID());
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
		    	
		    	ToolboxTool toolboxTool = new ToolboxTool(toolId, sadlDescription, parsedScript.SADL, parsedScript.code, parsedScript.source, resource, moduleDir.getName(), runtimeName);
		    	tools.put(toolId, toolboxTool);
		    	successfullyLoadedCount++;

	    		// Add to category, which gets sent to the client
	    		category.addTool(toolboxTool);
	    		
                // Set hidden if needed    		
                String hiddenStatus = "";
		    	if (toolHidden) {
		    		hiddenStatus = "HIDDEN";
		    		hiddenCount++;
		    	}

		    	logger.debug(String.format("loaded %s %s from %s %s" , toolId, sadlDescription.getName().getDisplayName(), resource, hiddenStatus));
		    }

		}

		// Update summary
		this.summary = "loaded " + moduleName + " " + successfullyLoadedCount + "/" + totalCount +
		" tools, " + disabledCount + " disabled, " + hiddenCount + " hidden";
		logger.info(summary);
	}

	public String getName() {
		return this.moduleName;
	}
	
	public String getNamePretty() {
		if (moduleName.equals("ngs")) {
			return moduleName.toUpperCase();
		} else {
			return StringUtils.capitalizeFirstLetter(moduleName);
		}
	}
	
	public List<ToolboxCategory> getCategories() {
		return this.categories;
	}
	
}