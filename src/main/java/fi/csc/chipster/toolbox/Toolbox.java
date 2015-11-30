package fi.csc.chipster.toolbox;

import java.io.File;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.chipster.toolbox.toolpartsparser.HeaderAsCommentParser;
import fi.csc.chipster.toolbox.toolpartsparser.JavaParser;
import fi.csc.chipster.toolbox.toolpartsparser.ToolPartsParser;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;


public class Toolbox {

	private static final Logger logger = Logger
			.getLogger(Toolbox.class);
	
	private List<ToolboxModule> modules = new LinkedList<ToolboxModule>();
	private File modulesDir;
	
	
	/**
	 * 
	 * @param the root workDir for the jobs of the computing service
	 * @throws IOException 
	 * @throws Exception
	 */
	public Toolbox(File modulesDir) throws IOException {
		this.modulesDir = modulesDir;
		loadModuleDescriptions();
	}
	
	public ToolboxTool getTool(String id) {
		
		// Iterate over modules and return description if it is found
		for (ToolboxModule module : modules) {
			ToolboxTool tool = module.getTool(id);
			if (tool != null) {
				return tool;
			}
		}
		
		// Nothing was found
		return null;
	}
	
	public List<ToolboxTool> getAll() {
		List<ToolboxTool> list = new LinkedList<ToolboxTool>();
		for (ToolboxModule module : modules) {
			list.addAll(module.getAll());
		}
		
		return list;
	}
	
	public List<ToolboxModule> getModules() {
		return this.modules;
	}
	
	
	/**
	 * @return a list of DescriptionMessages about available tool modules
	 * that can be sent to client.
	 */
	public List<ModuleDescriptionMessage> getModuleDescriptions() {
		
		LinkedList<ModuleDescriptionMessage> moduleDescriptions = new LinkedList<ModuleDescriptionMessage>();

		for (ToolboxModule module : modules) {
			moduleDescriptions.add(module.getModuleDescriptionMessage());
		}
		
	    return moduleDescriptions;
	}

	
	/**
	 * Load all the tool modules in this toolbox. Put them to the modules list.
	 * 
	 * @throws IOException
	 */
	private void loadModuleDescriptions() throws IOException {
		logger.info("loading modules");

		// Iterate over all module directories, and over all module files inside them
		List<String> moduleLoadSummaries = new LinkedList<String>();
		for (String moduleDirName : modulesDir.list()) {
			File moduleDir = new File(modulesDir, moduleDirName);

			if (moduleDir.isDirectory()) {

				// Load module specification files, if they exist (one module dir can contain multiple module specification files)
				for (String moduleFilename : moduleDir.list()) {
					if (moduleFilename.endsWith("-module.xml")) {
						File moduleFile = new File(moduleDir, moduleFilename);
						if (moduleFile.exists()) {

							// Load module
							logger.info("loading tools specifications from: " + moduleFilename);
							ToolboxModule module;
							String summary;
							try {
								module = new ToolboxModule(moduleDir, moduleFile);
								summary = module.getSummary();
							} catch (Exception e) {
								logger.warn("loading " + moduleFilename + " failed", e);
								continue;
							}
							// Register the module
							modules.add(module);
							moduleLoadSummaries.add(summary);
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
	}
	
	/**
	 * Toolbox modules use this to get the right parser for each runtime and tool type.
	 * 
	 * This is about separating the sadl part from for example R script, not parsing the sadl
	 * itself.
	 * 
	 * @param runtime
	 * @return
	 */
	static ToolPartsParser getToolPartsParser(String runtime) {
		if (runtime == null || runtime.isEmpty()) {
			return null;
		} else if (runtime.startsWith("python")) {
			return new HeaderAsCommentParser("#", runtime);
		} else if (runtime.startsWith("java")) {
			return new JavaParser();

		// add non-R stuff starting with R before this
		} else if (runtime.startsWith("R")) {
			return new HeaderAsCommentParser("#", runtime);
		} else {
			return null;
		}
	}
}