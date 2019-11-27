package fi.csc.chipster.toolbox;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.DirectoryStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.LinkedList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import fi.csc.chipster.toolbox.toolpartsparser.HeaderAsCommentParser;
import fi.csc.chipster.toolbox.toolpartsparser.JavaParser;
import fi.csc.chipster.toolbox.toolpartsparser.ToolPartsParser;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;

public class Toolbox {

	private static Logger logger = LogManager.getLogger();

	private List<ToolboxModule> modules = new LinkedList<ToolboxModule>();
	private byte[] zipContents;

	/**
	 * Loads tools.
	 * 
	 * @param modulesDir
	 * @param toolsBin
	 * @throws IOException
	 */
	public Toolbox(final Path modulesDir, File toolsBin) throws IOException {

		// load tools
		this.modules.addAll(loadModuleDescriptions(modulesDir, toolsBin));
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

	public ToolboxModule getModule(String name) {
		for (ToolboxModule module : modules) {
			if (module.getName().equals(name)) {
				return module;
			}
		}
		return null;
	}

	/**
	 * @return a list of DescriptionMessages about available tool modules that can
	 *         be sent to client.
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
	 * @param toolsBin2
	 * 
	 * @throws IOException
	 */
	public static List<ToolboxModule> loadModuleDescriptions(Path modulesDir, File toolsBin) throws IOException {

		// Iterate over all module directories, and over all module files inside
		// them
		List<ToolboxModule> modulesList = new LinkedList<ToolboxModule>();
		List<String> moduleLoadSummaries = new LinkedList<String>();
		int totalCount = 0;
		int successfullyLoadedCount = 0;
		int hiddenCount = 0;
		int disabledCount = 0;

		try (DirectoryStream<Path> modulesDirStream = Files.newDirectoryStream(modulesDir)) {
			for (Path moduleDir : modulesDirStream) {
				if (Files.isDirectory(moduleDir)) {

					// Load module specification files, if they exist (one
					// module dir can contain multiple module specification
					// files)
					try (DirectoryStream<Path> moduleDirStream = Files.newDirectoryStream(moduleDir, "*-module.xml")) {
						for (Path moduleFile : moduleDirStream) {

							// Load module
							logger.info("------ " + "loading tools specifications from: " + moduleFile + " ------");
							ToolboxModule module;
							String summary;
							try {
								module = new ToolboxModule(moduleDir, moduleFile, toolsBin);
								summary = module.getSummary();
							} catch (Exception e) {
								logger.warn("loading " + moduleFile + " failed", e);
								continue;
							}
							// Register the module
							modulesList.add(module);
							moduleLoadSummaries.add(summary);
							totalCount += module.getTotalCount();
							successfullyLoadedCount += module.getSuccessfullyLoadedCount();
							hiddenCount += module.getHiddenCount();
							disabledCount += module.getDisabledCount();

						}
					}
				}
			}
		}

		// print all summaries
		logger.info("------------ tool summary ------------");
		for (String summary : moduleLoadSummaries) {
			logger.info(summary);
		}
		logger.info("loaded total " + successfullyLoadedCount + "/" + totalCount + ", " + disabledCount + " disabled, "
				+ hiddenCount + " hidden");
		if (successfullyLoadedCount == totalCount) {
			logger.warn("****** ALL TOOLS LOADED SUCCESSFULLY ******");
		} else {
			logger.warn("****** LOADING SOME TOOLS FAILED ******");
		}

		logger.info("------------ tool summary ------------");

		return modulesList;
	}

	public void setZipContents(byte[] zipContents) {
		this.zipContents = zipContents;
	}

	public InputStream getZipStream() {
		return new ByteArrayInputStream(zipContents);
	}

	public byte[] getZipContents() {
		return zipContents;
	}

	/**
	 * Toolbox modules use this to get the right parser for each runtime and tool
	 * type.
	 * 
	 * This is about separating the sadl part from for example R script, not parsing
	 * the sadl itself.
	 * 
	 * @param toolId
	 * @return
	 */
	static ToolPartsParser getToolPartsParser(String toolId) {
		if (toolId == null || toolId.isEmpty()) {
			return null;
		}

		String defaultRuntime = getDefaultRuntime(toolId);
		if (defaultRuntime == null) {
			return null;
		}

		if (defaultRuntime.startsWith("python")) {
			return new HeaderAsCommentParser("#", RuntimeUtils.getToolDirFromRuntimeName(defaultRuntime));
		} else if (defaultRuntime.startsWith("java")) {
			return new JavaParser();

			// add non-R stuff starting with R before this
		} else if (defaultRuntime.startsWith("R")) {
			return new HeaderAsCommentParser("#", RuntimeUtils.getToolDirFromRuntimeName(defaultRuntime));
		} else {
			return null;
		}
	}

	static String getDefaultRuntime(String toolId) {
		if (toolId == null || toolId.isEmpty()) {
			return null;
		} else if (toolId.endsWith(".R")) {
			return "R";
		} else if (toolId.endsWith(".py")) {
			return "python";
		} else if (toolId.endsWith(".java")) {
			return "java";
		}
		return null;
	}
}