package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.microarray.util.XmlUtil;


/**
 * 
 * Any access to descriptions or visibleDescriptions should be use synchronized(this).
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
	private LinkedHashMap<String, AnalysisDescription> visibleDescriptions = new LinkedHashMap<String, AnalysisDescription>();
	private LinkedList<AnalysisHandler> handlers = new LinkedList<AnalysisHandler>(); 

	private File workDir;
	

	public ToolRepository(File workDir) {
		this.workDir = workDir;
	}
	
	
	public void addAnalysisHandler(AnalysisHandler handler) {
		handlers.add(handler);
	}
	
	
	public void loadOperation(String sourceResourceName, boolean hidden) throws AnalysisException {
		AnalysisDescription description = loadDescription(sourceResourceName);
		if (description != null) {
			addDescription(description, hidden);
		}
	}
	
	
	private AnalysisDescription loadDescription(String sourceResourceName) throws AnalysisException {
		for (AnalysisHandler handler : handlers) {
			logger.debug("using handler " + handler.getClass().getSimpleName() + ", checking " + sourceResourceName);
			if (handler.supports(sourceResourceName)) {
				return handler.handle(sourceResourceName);
			}
		}
		logger.warn("none of the loaded handlers support " + sourceResourceName);
		return null;
	}
	
	
	
	private void addDescription(AnalysisDescription description, boolean hidden) {
		logger.debug("added operation " + description.getFullName());
		
		synchronized(this) {
			descriptions.put(description.getFullName(), description);
			if (!hidden) {
				visibleDescriptions.put(description.getFullName(), description);
			}
		}
	}
	
	public AnalysisDescription getDescription(String fullName) throws AnalysisException {
		AnalysisDescription desc; 

		// get the description
		synchronized(this) {
			desc = descriptions.get(fullName);
		}
		
		// check if description needs to be updated
		if (desc != null && !desc.isUptodate()) {
			AnalysisDescription newDescription = loadDescription(desc.getSourceResourceName());
			if (newDescription != null) {
				synchronized(this) {
					descriptions.remove(fullName);
					descriptions.put(newDescription.getFullName(), newDescription);
					assert(newDescription.getFullName().equals(fullName));
					if (visibleDescriptions.containsKey(fullName)) {
						visibleDescriptions.remove(fullName);
						visibleDescriptions.put(newDescription.getFullName(), newDescription);
					}
				}
				desc = newDescription;
			}
		}

		return desc; 
	}
	
	/**
	 * Returns one huge VVSADL block that contains all loaded analysis 
	 * descriptions.
	 * @return huge block
	 */
	public StringBuffer serialiseAsStringBuffer() {
		StringBuffer buf = new StringBuffer();
		for (AnalysisDescription description : visibleDescriptions.values()) {
			buf.append(description.getVVSADL());
		}
		return buf;
	}

	private void loadToolPackages() throws FileNotFoundException, SAXException, IOException, ParserConfigurationException, IllegalArgumentException, SecurityException, InstantiationException, IllegalAccessException, InvocationTargetException, NoSuchMethodException, ClassNotFoundException { 
		logger.info("loading tool packages");
		
		// FIXME add real path
		File toolPackageConfig = new File("tool-packages.xml");
		
		Document document = XmlUtil.getInstance().parseReader(new FileReader(toolPackageConfig));
		Element toolPackagesElement = (Element)document.getElementsByTagName("runtimes").item(0);

		for (Element toolPackageElement: XmlUtil.getChildElements(toolPackagesElement, "runtime")) {
			logger.info("loading tool package: " + XmlUtil.getChildElement(toolPackageElement, "name").getTextContent());
			Element handlerElement = XmlUtil.getChildElement(toolPackageElement, "handler");
			
			String handlerClassName = XmlUtil.getChildElement(handlerElement, "class").getTextContent();

			// parameters to handler
			HashMap<String, String> parameters = new HashMap<String, String>();

			// comp work dir
			parameters.put("workDir", this.workDir.toString());

			// parameters from tool packages config
			for (Element parameterElement: XmlUtil.getChildElements(handlerElement, "parameter")) {
				String paramName = XmlUtil.getChildElement(parameterElement, "name").getTextContent();
				String paramValue = XmlUtil.getChildElement(parameterElement, "value").getTextContent(); 
				parameters.put(paramName, paramValue);
			}
			
			// instantiate handler
			AnalysisHandler handler = (AnalysisHandler)Class.forName(handlerClassName).getConstructor(HashMap.class).newInstance(parameters);
			// FIXME maybe not needed at all?
			//descriptionRepository.addAnalysisHandler(handler);

			// load tools
			Element toolsElement = XmlUtil.getChildElement(toolPackageElement, "tools"); 
			for (Element toolElement: XmlUtil.getChildElements(toolsElement, "tool")) {
			
			}
		}
		
//		for (int i = 0; i < toolPackages.getLength(); i++) {
//			Element toolPackageElement = (Element)toolPackages.item(i);
//			System.out.println(toolPackageElement.getElementsByTagName("name").item(0).getTextContent());
//		}
	}

	
}
