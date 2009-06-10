package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.LinkedHashMap;

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
 * description to be updated and being briefly unavailable during the
 * update.
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
			AnalysisDescription newDescription = desc.getHandler().handle(desc.getSourceResourceName());
			if (newDescription != null) {

				// name (id) of the tool has not changed
				if (desc.getFullName().equals(newDescription.getFullName())) {
					descriptions.remove(fullName);
					descriptions.put(newDescription.getFullName(), newDescription);
					if (supportedDescriptions.containsKey(desc.getFullName())) {
						supportedDescriptions.remove(desc.getFullName());
						supportedDescriptions.put(newDescription.getFullName(), newDescription);
					}
					if (visibleDescriptions.containsKey(desc.getFullName())) {
						visibleDescriptions.remove(desc.getFullName());
						visibleDescriptions.put(newDescription.getFullName(), newDescription);
					}
					return newDescription;
				} 

				// name (id) of the tool has changed
				else {
					logger.warn("name of the tool has changed after loading from custom-scripts, keeping both old and new");
					if (descriptions.containsKey(newDescription.getFullName())){
						logger.warn("descriptions already contains a tool with the new name, ignoring custom-scripts");
						return desc;
					} 
					// add the tool with the new name
					descriptions.put(newDescription.getFullName(), newDescription);
					if (supportedDescriptions.containsKey(desc.getFullName())) {
						supportedDescriptions.put(newDescription.getFullName(), newDescription);
					}
					if (visibleDescriptions.containsKey(desc.getFullName())) {
						visibleDescriptions.put(newDescription.getFullName(), newDescription);
					}
					return newDescription;
				}
			}
		}
		return desc; 
	}
	
	/**
	 * Returns one huge VVSADL block that contains all loaded analysis 
	 * descriptions.
	 * @return huge block
	 */
	public synchronized StringBuffer serialiseAsStringBuffer() {
		StringBuffer buf = new StringBuffer();
		for (AnalysisDescription description : visibleDescriptions.values()) {
			buf.append(description.getVVSADL());
		}
		return buf;
	}

	public synchronized boolean supports(String fullName) {
		return supportedDescriptions.containsKey(fullName);
	}

	private void loadRuntimes(File workDir) throws IllegalArgumentException, SecurityException, InstantiationException, IllegalAccessException, InvocationTargetException, NoSuchMethodException, ClassNotFoundException, IOException, SAXException, ParserConfigurationException  { 
		logger.info("loading runtimes");

		File runtimeConfig = new File(DirectoryLayout.getInstance().getConfDir(), "runtimes.xml");

		Document document = XmlUtil.getInstance().parseReader(new FileReader(runtimeConfig));
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
				logger.info("runtime " + runtimeName + " disabled since handler is disabled");
			}

			// add to runtimes
			ToolRuntime runtime = new ToolRuntime(runtimeName, handler, runtimeDisabled); 
			this.runtimes.put(runtimeName, runtime);
		}
	}	

	private void loadTools() throws IOException, SAXException, ParserConfigurationException { 
		logger.info("loading tools");
		
		File toolConfig = new File(DirectoryLayout.getInstance().getConfDir(), "tools.xml");
		
		Document document = XmlUtil.getInstance().parseReader(new FileReader(toolConfig));
		Element toolsElement = (Element)document.getElementsByTagName("tools").item(0);

		
		for (Element toolElement: XmlUtil.getChildElements(toolsElement, "tool")) {

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
				logger.warn("could not find runtime " + runtimeName + " for " + sourceResourceName);
				continue;
			}
			
			AnalysisDescription description;
			try {
				description = runtime.getHandler().handle(sourceResourceName);
			} catch (AnalysisException e) {
				logger.warn("could not create description for " + sourceResourceName);
				continue;
			}

			// add to descriptions
			if (descriptions.containsKey(description.getFullName())) {
				logger.warn("description with the name " + description.getFullName() + " already exists, keeping the original");
				continue;
			}
			descriptions.put(description.getFullName(), description);
			String disabledStatus = "";
			if (!runtime.isDisabled() && !toolDisabled) {
				supportedDescriptions.put(description.getFullName(), description);
			} else {
				disabledStatus = " DISABLED";
			}
			String hiddenStatus = "";
			if (!toolHidden) {
				visibleDescriptions.put(description.getFullName(), description);
			} else {
				hiddenStatus = " HIDDEN";
			}
			
			logger.info("loaded " + description.getFullName().replace("\"", "") + " " + description.getSourceResourceFullPath() + disabledStatus + hiddenStatus);
		}
	}
}
