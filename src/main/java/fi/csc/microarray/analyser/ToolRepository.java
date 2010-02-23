package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.config.DirectoryLayout;
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
 * @author hupponen
 *
 */
public class ToolRepository {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(ToolRepository.class);
	
	private LinkedHashMap<String, AnalysisDescription> descriptions = new LinkedHashMap<String, AnalysisDescription>(); 
	private LinkedHashMap<String, AnalysisDescription> supportedDescriptions = new LinkedHashMap<String, AnalysisDescription>();
	private LinkedHashMap<String, AnalysisDescription> visibleDescriptions = new LinkedHashMap<String, AnalysisDescription>();
	
	private HashMap<String, ToolRuntime> runtimes = new HashMap<String, ToolRuntime>();
	
	
	/**
	 * 
	 * @param the root workDir for the jobs of the computing service
	 * @throws Exception
	 */
	public ToolRepository(File workDir) throws Exception {
		loadRuntimes(workDir);
		loadTools();
	}
	
	
	
	
	public synchronized AnalysisDescription getDescription(String fullName) throws AnalysisException {
		AnalysisDescription desc; 

		// get the description
		desc = descriptions.get(fullName);

		// check if description needs to be updated
		if (desc != null && !desc.isUptodate()) {
			updateDescription(desc);
		}
		
		// return the possibly updated description
		return descriptions.get(fullName); 
	}




	
	/**
	 * Returns one huge VVSADL block that contains all loaded analysis 
	 * descriptions.
	 * @return huge block
	 * @throws AnalysisException 
	 */
	public synchronized StringBuffer serialiseAsStringBuffer() throws AnalysisException {
		StringBuffer buf = new StringBuffer();

		// find descs that need to be updated (custom script available)
		List<AnalysisDescription> descsToBeUpdated = new LinkedList<AnalysisDescription>();
		for (AnalysisDescription description : visibleDescriptions.values()) {
			if (!description.isUptodate()) {
				descsToBeUpdated.add(description);
			}
		}
		
		// update (can't update in the previous loop, would cause concurrent modification)
		for (AnalysisDescription description: descsToBeUpdated) {
			updateDescription(description);
		}
		
		// get the descriptions
		for (AnalysisDescription description: visibleDescriptions.values()) {
			buf.append(description.getVVSADL());
		}
		return buf;
	}

	public synchronized boolean supports(String fullName) {
		return supportedDescriptions.containsKey(fullName);
	}

	
	private void updateDescription(AnalysisDescription desc) throws AnalysisException {
		AnalysisDescription newDescription = desc.getHandler().handle(desc.getSourceResourceName());
		if (newDescription != null) {
			newDescription.setUpdatedSinceStartup();
			
			// name (id) of the tool has not changed
			if (desc.getFullName().equals(newDescription.getFullName())) {
				
				// replace the old description with the same name
				descriptions.put(newDescription.getFullName(), newDescription);
				if (supportedDescriptions.containsKey(desc.getFullName())) {
					supportedDescriptions.put(newDescription.getFullName(), newDescription);
				}
				if (visibleDescriptions.containsKey(desc.getFullName())) {
					visibleDescriptions.put(newDescription.getFullName(), newDescription);
				}
			} 

			// name (id) of the tool has changed
			else {
				logger.warn("name of the tool has changed after loading from custom-scripts, keeping both old and new");
				if (descriptions.containsKey(newDescription.getFullName())){
					logger.warn("descriptions already contains a tool with the new name, ignoring custom-scripts");
					return;
				} 
				// add the tool with the new name
				descriptions.put(newDescription.getFullName(), newDescription);
				if (supportedDescriptions.containsKey(desc.getFullName())) {
					supportedDescriptions.put(newDescription.getFullName(), newDescription);
				}
				if (visibleDescriptions.containsKey(desc.getFullName())) {
					visibleDescriptions.put(newDescription.getFullName(), newDescription);
				}
			}
		}
	}

	
	
	private void loadRuntimes(File workDir) throws IllegalArgumentException, SecurityException, InstantiationException, IllegalAccessException, InvocationTargetException, NoSuchMethodException, ClassNotFoundException, IOException, SAXException, ParserConfigurationException  { 
		logger.info("loading runtimes");

		File runtimeConfig = new File(DirectoryLayout.getInstance().getConfDir(), "runtimes.xml");

		Document document = XmlUtil.parseReader(new FileReader(runtimeConfig));
		Element runtimesElement = (Element)document.getElementsByTagName("runtimes").item(0);


		for (Element runtimeElement: XmlUtil.getChildElements(runtimesElement, "runtime")) {
			String runtimeName = XmlUtil.getChildElement(runtimeElement, "name").getTextContent();
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
			String handlerClassName = XmlUtil.getChildElement(handlerElement, "class").getTextContent();

			// parameters to handler
			HashMap<String, String> parameters = new HashMap<String, String>();

			// comp work dir
			parameters.put("workDir", workDir.toString());

			// parameters from config
			for (Element parameterElement: XmlUtil.getChildElements(handlerElement, "parameter")) {
				String paramName = XmlUtil.getChildElement(parameterElement, "name").getTextContent();
				String paramValue = XmlUtil.getChildElement(parameterElement, "value").getTextContent(); 
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

	private void loadTools() throws IOException, SAXException, ParserConfigurationException { 
		logger.info("loading tools");
		
		String [] toolFiles = new String[] { "tools.xml", "emboss-tools.xml" };
		for (String toolFileName: toolFiles) {
			File toolFile = new File(DirectoryLayout.getInstance().getConfDir(), toolFileName);
			if (toolFile.exists()) {
				logger.info("loading from " + toolFileName);
				loadToolsFromFile(toolFile);
			}
		}
	}


	private void loadToolsFromFile(File toolFile) throws FileNotFoundException, SAXException, IOException, ParserConfigurationException {
		File toolConfig = toolFile;

		Document document = XmlUtil.parseReader(new FileReader(toolConfig));
		Element toolsElement = (Element)document.getElementsByTagName("tools").item(0);

		int totalCount = 0;
		int successfullyLoadedCount = 0;
		int hiddenCount = 0;
		int disabledCount = 0;
		for (Element toolElement: XmlUtil.getChildElements(toolsElement, "tool")) {
			totalCount++;

			// tool name
			String sourceResourceName = toolElement.getTextContent();
			logger.debug("loading " + sourceResourceName);

			// runtime
			String runtimeName = toolElement.getAttribute("runtime");

			// disabled or hidden?
			boolean toolDisabled = toolElement.getAttribute("disabled").equals("true");
			boolean toolHidden = toolElement.getAttribute("hidden").equals("true");

			// load the tool
			ToolRuntime runtime = runtimes.get(runtimeName);
			if (runtime == null) {
				logger.warn("loading " + sourceResourceName + " failed, could not find runtime " + runtimeName);
				continue;
			}

			AnalysisDescription description;
			try {
				description = runtime.getHandler().handle(sourceResourceName);
			} catch (Exception e) {
				logger.warn("loading " + sourceResourceName + " failed, could not create description", e);
				continue;
			}

			// add to descriptions
			if (descriptions.containsKey(description.getFullName())) {
				logger.warn("loading " + sourceResourceName + " failed, description with the name " + description.getFullName() + " already exists");
				continue;
			}
			descriptions.put(description.getFullName(), description);
			successfullyLoadedCount++;

			String disabledStatus = "";
			if (!runtime.isDisabled() && !toolDisabled) {
				supportedDescriptions.put(description.getFullName(), description);
			} else {
				disabledStatus = " DISABLED";
				disabledCount++;
			}
			String hiddenStatus = "";
			if (!toolHidden) {
				visibleDescriptions.put(description.getFullName(), description);
			} else {
				hiddenStatus = " HIDDEN";
				hiddenCount++;
			}

			logger.info("loaded " + description.getFullName().replace("\"", "") + " " + description.getSourceResourceFullPath() + disabledStatus + hiddenStatus);
		}
		logger.info("loaded " + successfullyLoadedCount + "/" + totalCount + " tools, " + disabledCount + " disabled, " + hiddenCount + " hidden");

	}

}