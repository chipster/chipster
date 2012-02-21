package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
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
import fi.csc.microarray.util.XmlUtil;


/**
 * <p>Manages analysis tools, modules and runtimes. Analysis tools
 * are visible to AnalyserServer as ToolDescriptions, whereas runtimes and modules
 * are not visible outside of this class.</p>
 *
 * <p>Access is strictly synchronised, because all operations
 * may lead to module or script updates if files on the disk have changed. To avoid
 * deadlocking, dependencies must be kept one way: RepositoryModule never calls ToolRepository and
 * ToolDescription never calls RepositoryModule.</p>
 *
 * <p>Initialisation of the repository can fail, but after that it is robust. The repository can be
 * reload while it is running, but if reloading fails, it retains the original state. In other
 * words, only constructor throws exceptions. Other errors are just logged.</p>
 *   
 * @see RepositoryModule
 * @see ToolDescription
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
	private List<RepositoryModule> modules = new LinkedList<RepositoryModule>();
		
	/**
	 * 
	 * @param the root workDir for the jobs of the computing service
	 * @throws Exception
	 */
	public ToolRepository(File workDir) throws AnalysisException {
		loadRuntimes(workDir);
		loadModuleDescriptions();
	}
	
	public synchronized ToolDescription getDescription(String id) {
		
		// Iterate over modules and return description if it is found
		for (RepositoryModule module : modules) {
			ToolDescription desc = module.getDescription(id);
			
			if (desc != null) {
				return desc;
			}
		}
		
		// Nothing was found
		return null;
	}
	
	/**
	 * @return true if this comp service can run the given tool.
	 */
	public synchronized boolean supports(String toolId) {
		for (RepositoryModule module : modules) {
			if (module.isSupportedDescription(toolId)) {
				return true;
			}
		}
		return false;
	}

	
	/**
	 * @return a list of DescriptionMessages about available modules
	 * that can be sent to client.
	 * @throws ParserConfigurationException
	 * @throws SAXException
     * @throws IOException
	 */
	public synchronized List<ModuleDescriptionMessage> getModuleDescriptions() {
		
		LinkedList<ModuleDescriptionMessage> moduleDescriptions = new LinkedList<ModuleDescriptionMessage>();

		for (RepositoryModule module : modules) {
			moduleDescriptions.add(module.getModuleDescriptionMessage());
		}
		
	    return moduleDescriptions;
	}

	
	/**
	 * Load available runtimes.
	 * 
	 * @param workDir
	 * @throws AnalysisException 
	 */
	private synchronized void loadRuntimes(File workDir) throws AnalysisException  { 
		logger.info("loading runtimes");

		try {
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

		} catch (Exception e) {
			throw new AnalysisException(e);
		}
	}

	/**
	 * Generate a list containing information about all modules
	 * available in this analyser server.
	 * 
	 * @throws IOException
	 * @throws SAXException
	 * @throws ParserConfigurationException
	 */
	private synchronized void loadModuleDescriptions()
	       throws AnalysisException {
		logger.info("loading modules");

		try {
			// Iterate over all module directories, and over all module files inside them
			List<String> moduleLoadSummaries = new LinkedList<String>();
			for (String moduleDirName : DirectoryLayout.getInstance().getModulesDir().list()) {
				File moduleDir = new File(DirectoryLayout.getInstance().getModulesDir(), moduleDirName);

				if (moduleDir.isDirectory()) {

					// Load module specification files, if they exist (one module dir can contain multiple module specification files)
					for (String moduleFilename : moduleDir.list()) {
						if (moduleFilename.endsWith("-module.xml")) {
							File moduleFile = new File(moduleDir, moduleFilename);
							if (moduleFile.exists()) {

								// Load module
								logger.info("loading tools specifications from: " + moduleFilename);
								RepositoryModule module = new RepositoryModule(moduleDir, moduleFile, runtimes);

								// Register the module
								modules.add(module);
								moduleLoadSummaries.add(module.getSummary());
							}
						}
					}
				}
			}

			// print all summaries
			logger.info("------ tool summary ------ ");
			for (String summary : moduleLoadSummaries) {
				logger.info(summary);
			}
			logger.info("------ tool summary ------ ");

		} catch (Exception e) {
			throw new AnalysisException(e);
		}
	}
}