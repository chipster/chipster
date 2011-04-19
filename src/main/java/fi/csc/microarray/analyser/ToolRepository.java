package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage.Category;
import fi.csc.microarray.util.XmlUtil;


/**
 * The ToolRepository manages analysis tools and runtimes. Analysis tools
 * are visible to AnalyserServer as AnalysisDescriptions, where as runtimes
 * are not visible outside of this class.
 * 
 * After initialization, any access to descriptions maps should be use 
 * synchronized using this, since reading descriptions may cause a
 * description to be updated.
 *  
 * @author Taavi Hupponen, Aleksi Kallio
 *
 */
public class ToolRepository {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(ToolRepository.class);
	
	private HashMap<String, ToolRuntime> runtimes = new HashMap<String, ToolRuntime>();
	private List<ModuleDescriptionMessage> moduleDescriptions = new LinkedList<ModuleDescriptionMessage>();
	private List<RepositoryModule> modules = new LinkedList<RepositoryModule>();
		
	/**
	 * 
	 * @param the root workDir for the jobs of the computing service
	 * @throws Exception
	 */
	public ToolRepository(File workDir) throws Exception {
		loadRuntimes(workDir);
		loadModuleDescriptions();
	}
	
	public synchronized AnalysisDescription getDescription(String id) throws AnalysisException {
		AnalysisDescription desc = null; 
		RepositoryModule moduleWithDesc = null;
		
		// Get the description
		for (RepositoryModule module : modules) {
			if (module.hasDescription(id)) {
				desc = module.getDescription(id);
				moduleWithDesc = module;
				break;
			}
		}
		
		// Return null if nothing is found
		if (desc == null) {
			return null;
		}

		// Check if description needs to be updated
		if (desc != null && !desc.isUptodate()) {
			moduleWithDesc.updateDescription(desc);
		}
		
		// Return the possibly updated description
		return moduleWithDesc.getDescription(id); 
	}
	
	public synchronized boolean supports(String id) {
		for (RepositoryModule module : modules) {
			if (module.isSupportedDescription(id)) {
				return true;
			}
		}
		return false;
	}
	

	/**
	 * Load available runtimes.
	 * 
	 * @param workDir
	 */
	private void loadRuntimes(File workDir) throws IllegalArgumentException, SecurityException, InstantiationException, IllegalAccessException, InvocationTargetException, NoSuchMethodException, ClassNotFoundException, IOException, SAXException, ParserConfigurationException  { 
		logger.info("loading runtimes");

		File runtimeConfig = new File(DirectoryLayout.getInstance().getConfDir(), "runtimes.xml");

		Document document = XmlUtil.parseReader(new FileReader(runtimeConfig));
		Element runtimesElement = (Element)document.getElementsByTagName("runtimes").item(0);


		for (Element runtimeElement: XmlUtil.getChildElements(runtimesElement, "runtime")) {
			String runtimeName = XmlUtil.getChildElement(runtimeElement, "name").getTextContent().trim();
			logger.info("loading runtime " + runtimeName);
			if (runtimes.containsKey(runtimeName)) {
				logger.warn("runtime with the same name " + runtimeName + " already loaded, keeping the first one");
				continue;
			}
			
			
			boolean runtimeDisabled = runtimeElement.getAttribute("disabled").equals("true");
			if (runtimeDisabled) {
				logger.info("runtime " + runtimeName + " disabled by config");
			}
			
			// handler
			Element handlerElement = XmlUtil.getChildElement(runtimeElement, "handler");
			String handlerClassName = XmlUtil.getChildElement(handlerElement, "class").getTextContent().trim();

			// parameters to handler
			HashMap<String, String> parameters = new HashMap<String, String>();

			// comp work dir
			parameters.put("workDir", workDir.toString());

			// parameters from config
			for (Element parameterElement: XmlUtil.getChildElements(handlerElement, "parameter")) {
				String paramName = XmlUtil.getChildElement(parameterElement, "name").getTextContent().trim();
				String paramValue = XmlUtil.getChildElement(parameterElement, "value").getTextContent().trim(); 
				parameters.put(paramName, paramValue);
			}

			// instantiate handler
			AnalysisHandler handler = (AnalysisHandler)Class.forName(handlerClassName).getConstructor(HashMap.class).newInstance(parameters);
			
			// disabled
			if (handler.isDisabled()) {
				runtimeDisabled = true;
				logger.info("runtime " + runtimeName + " disabled as handler is disabled");
			}

			// add to runtimes
			ToolRuntime runtime = new ToolRuntime(runtimeName, handler, runtimeDisabled); 
			this.runtimes.put(runtimeName, runtime);
		}
	}
	
	/**
	 * @return a list of DescriptionMessages about available modules
	 * that can be sent to client.
	 * @throws ParserConfigurationException
	 * @throws SAXException
     * @throws IOException
	 */
	public List<ModuleDescriptionMessage> getModuleDescriptions()
	        throws ParserConfigurationException, SAXException, IOException {
	    return moduleDescriptions;
	}

	/**
	 * Generate a list containing information about all modules
	 * available in this analyser server.
	 * 
	 * @throws IOException
	 * @throws SAXException
	 * @throws ParserConfigurationException
	 */
	private void loadModuleDescriptions()
	       throws IOException, SAXException, ParserConfigurationException {
		logger.info("loading modules");

		// Iterate over all module directories, and over all module files inside them
		for (String moduleDirName : DirectoryLayout.getInstance().getModulesDir().list()) {
			File moduleDir = new File(DirectoryLayout.getInstance().getModulesDir(), moduleDirName);
			if (moduleDir.isDirectory()) {
				
				
				// Make directory as a known source of tools
				// FIXME
				logger.info("detected tool directory: " + moduleDir);
				
				// Load module specification files, if they exist
				for (String moduleFilename : moduleDir.list()) {
				    if (moduleFilename.endsWith("-module.xml")) {
			            File moduleFile = new File(moduleDir, moduleFilename);
			            if (moduleFile.exists()) {
			                logger.info("loading tools specifications from: " + moduleFilename);
			                loadModule(moduleDir, moduleFile);
			            }
				    }
				}
			}
		}
	}

	/**
	 * Parses a module file and loads all tools listed in it.
	 * 
	 * During parsing also construct description message that
	 * can be sent to the client.
	 * 
	 * Disabled tools are not available in this analyser instance,
	 * but are available in some other instances running, so the
	 * client still has to be informed about them.
	 * 
	 * Hidden tools are available in this instance, but are not
	 * shown in the normal operation list in user interface. Still,
	 * we should pass them to the client.
	 * 
	 * @param toolFile
	 * @throws FileNotFoundException
	 * @throws SAXException
	 * @throws IOException
	 * @throws ParserConfigurationException
	 */
	private void loadModule(File moduleDir, File toolFile)
	    throws FileNotFoundException, SAXException,
	           IOException, ParserConfigurationException {

		Document document = XmlUtil.parseReader(new FileReader(toolFile));
		Element moduleElement = (Element)document.getElementsByTagName("module").item(0);
		
		// module name 
		String moduleName = moduleElement.getAttribute("name");
		if (moduleName.isEmpty()) {
			logger.warn("not loading a module without a name");
			return;
		}
		
		// construct the module
		RepositoryModule module = new RepositoryModule(moduleDir, toolFile, moduleName);
		
		// stats
	    int totalCount = 0;
		int successfullyLoadedCount = 0;
	    int hiddenCount = 0;
	    int disabledCount = 0;
		
		// load categories
		for (Element categoryElement: XmlUtil.getChildElements(moduleElement, "category")) {

			// name
			String categoryName = categoryElement.getAttribute("name");
			if (categoryName.isEmpty()) {
				logger.warn("not loading a category without a name");
				continue;
			}
			
			// enabled or disabled
			boolean categoryDisabled = categoryElement.getAttribute("disabled").equals("true");
			if (categoryDisabled) {
				logger.info("not loading category " + categoryName + ": disabled");
				continue;
			}
			
			// color
			String categoryColor = categoryElement.getAttribute("color");
			if (categoryColor.isEmpty()) {
				logger.warn("not loading category " + categoryName + ": no color");
				continue;
			}
			
		    // visibility
		    boolean categoryHidden = Boolean.valueOf(categoryElement.getAttribute("hidden"));

		    // add category to the module description message
		    Category category = new Category(categoryName, categoryColor, categoryHidden);
		    
		    // load tools
		    for (Element toolElement: XmlUtil.getChildElements(categoryElement, "tool")) {
		    	totalCount++;

		    	Element resourceElement = XmlUtil.getChildElement(toolElement, "resource");
		    	if (resourceElement == null) {
		    		logger.warn("not loading a tool without resource element");
		    		continue;
		    	}
		    	String toolFilename = resourceElement.getTextContent().trim();
		    	if (toolFilename == null || toolFilename.isEmpty()) {
		    		logger.warn("not loading a tool with empty filename");
		    		continue;
		    	}

		    	// tool visibility
		    	// currently tools can be hidden only if their category is hidden
		    	boolean toolDisabled = toolElement.getAttribute("disabled").equals("true");
		    	boolean toolHidden = categoryHidden;

		    	// runtime
		    	String runtimeName = toolElement.getAttribute("runtime");
		    	ToolRuntime runtime = runtimes.get(runtimeName);
		    	if (runtime == null) {
		    		logger.warn("not loading " + toolFilename + ": runtime " + runtimeName + " not found");
		    		continue;
		    	}

		    	// parameters
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

		    		// parameter ok
		    		parameters.put(parameterName, parameterValue);
		    	}
		    	if (!parametersOk) {
		    		logger.warn("not loading " + toolFilename + ": parameters not ok");
		    		continue;
		    	}
		    	
		    	// create the analysis description
		    	AnalysisDescription description;
		    	try {
		    	    // checked cached descriptions
		    	    // cached descriptions are updated when needed
		    	    AnalysisDescription cachedDescription = getDescription(toolFilename);
		    	    if (cachedDescription != null) {
		    	        description = cachedDescription;
		    	    } else {
		                description = runtime.getHandler().handle(moduleDir, toolFilename, parameters);
		    	    }
		    	    
		    	} catch (Exception e) {
		    		logger.warn("loading " + toolFilename + " failed, could not create description", e);
		    		continue;
		    	}
		    	module.putDescription(description.getID(), description);

		    	successfullyLoadedCount++;

		    	// disabled or not
		    	String disabledStatus = "";
		    	if (!runtime.isDisabled() && !toolDisabled) {
		    		// add to supported descriptions list
		    		module.markSupportedDescription(description.getID());
		    	} else {
		    		disabledStatus = " DISABLED";
		    		disabledCount++;
		    	}

	    		// add to category, which gets sent to the client
	    		category.addTool(description.getSADL(), description.getHelpURL());
	    		
                // hidden or not    		
                String hiddenStatus = "";
		    	if (toolHidden) {
		    		hiddenStatus = " HIDDEN";
		    		hiddenCount++;
		    	}

		    	logger.info("loaded " + description.getID() + " " + description.getDisplayName() + " " +
		    			description.getSourceResourceFullPath() + disabledStatus + hiddenStatus);
		    }

		    // add the category to the module description message
		    module.addDescriptionCategory(category);
		}

		logger.info("loaded " + moduleName + " " + successfullyLoadedCount + "/" + totalCount +
				" tools, " + disabledCount + " disabled, " + hiddenCount + " hidden");

		// add to modules
		modules.add(module);
		moduleDescriptions.add(module.getModuleDescriptionMessage());
	}
}