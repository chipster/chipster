package fi.csc.microarray.comp;

import java.io.File;
import java.io.FileReader;
import java.util.HashMap;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.XmlUtil;


public class RuntimeRepository {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(RuntimeRepository.class);
	
	private HashMap<String, ToolRuntime> runtimes = new HashMap<String, ToolRuntime>();
		
	/**
	 * 
	 * @param the root workDir for the jobs of the computing service
	 * @throws Exception
	 */
	public RuntimeRepository(File workDir) throws CompException {
		loadRuntimes(workDir);
	}
	
	public ToolRuntime getRuntime(String name) {
		return runtimes.get(name);
	}
	
	
	/**
	 * Load available runtimes.
	 * 
	 * @param workDir
	 * @throws CompException 
	 */
	private synchronized void loadRuntimes(File workDir) throws CompException  { 
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

				// instantiate job factory
				JobFactory jobFactory = (JobFactory)Class.forName(handlerClassName).getConstructor(HashMap.class).newInstance(parameters);

				// disabled
				if (jobFactory.isDisabled()) {
					runtimeDisabled = true;
					logger.info("runtime " + runtimeName + " disabled as handler is disabled");
				}

				// add to runtimes
				ToolRuntime runtime = new ToolRuntime(runtimeName, jobFactory, runtimeDisabled); 
				this.runtimes.put(runtimeName, runtime);
			}

		} catch (Exception e) {
			throw new CompException(e);
		}
	}
}