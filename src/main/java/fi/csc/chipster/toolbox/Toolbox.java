package fi.csc.chipster.toolbox;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import de.schlichtherle.truezip.file.TFile;
import fi.csc.chipster.toolbox.toolpartsparser.HeaderAsCommentParser;
import fi.csc.chipster.toolbox.toolpartsparser.JavaParser;
import fi.csc.chipster.toolbox.toolpartsparser.ToolPartsParser;
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.util.Files;


public class Toolbox {

	private static final Logger logger = Logger
			.getLogger(Toolbox.class);
	
	private static final String MODULES_DIR_NAME = "modules";
	private static final String TOOLS_DIST_BASENAME = "chipster-tools";
	
	private List<ToolboxModule> modules = new LinkedList<ToolboxModule>();
	private File modulesDir = null;
	
	
	/**
	 * 
	 * @param the root workDir for the jobs of the computing service
	 * @throws IOException 
	 * @throws URISyntaxException 
	 * @throws Exception
	 */
	public Toolbox(final File modulesDir) throws IOException, URISyntaxException {
		// decide tools location
		File rootDir = new TFile(".").getCanonicalFile();
		
		findModulesDir(modulesDir, rootDir);
		
		
		// load tools
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
	
	public ToolboxModule getModule(String name) {
		for (ToolboxModule module : modules) {
			if (module.getName().equals(name)) {
				return module;
			}
		}
		return null;
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
	 * Try given modulesDir first. 
	 * 
	 * If no luck, search rootDir for chipster-tools dist, could be dir or tar.gz.
	 * 
	 * @param modulesDir
	 * @param rootDir
	 * @throws IOException
	 * @throws FileNotFoundException
	 * @throws URISyntaxException 
	 */
	private void findModulesDir(final File modulesDir, File rootDir) throws IOException, FileNotFoundException, URISyntaxException {
		// modules dir
		logger.info("looking for modules dir " + modulesDir.getAbsolutePath());
		if (modulesDir.exists() && modulesDir.isDirectory()) {
			this.modulesDir = modulesDir;
			logger.info("modules dir " + modulesDir.getAbsolutePath() + " found");
		}
		
		// chipster-tools-x.y.z dir
		else {
			logger.info("modules dir " + modulesDir.getAbsolutePath() + " not found");
			logger.info("looking for " + TOOLS_DIST_BASENAME + "-x.y.z dir in " + rootDir);
			
			File toolsDir = Files.getLatestVersion(rootDir, TOOLS_DIST_BASENAME, null);
			if (toolsDir != null && toolsDir.exists()) {
				logger.info("found " + toolsDir);
				logger.info("looking for " + toolsDir.getAbsolutePath() + File.separator + MODULES_DIR_NAME);
				File possibleModulesDir = new File(toolsDir, MODULES_DIR_NAME); 
				if (possibleModulesDir.exists() && possibleModulesDir.isDirectory()) {
					this.modulesDir = possibleModulesDir;
					logger.info("modules dir " + this.modulesDir + " found");
				} else {
					logger.info("modules dir " + possibleModulesDir + " not found");
				}
			} else {
				logger.info(rootDir + File.separator + TOOLS_DIST_BASENAME + "-x.y.z dir not found");
			}
			
			// still not found, try chipster-tools-x.y.z-tar.gz
			if (this.modulesDir == null) {
				logger.info("looking for " + new File(rootDir, TOOLS_DIST_BASENAME + "-x.y.z" + ".tar.gz"));
				File toolsTar = Files.getLatestVersion(rootDir, TOOLS_DIST_BASENAME, "tar.gz");
				loadFromTar(toolsTar);
			}

			// still not found, try chipster-tools-x.y.z-tar.gz from classpath
			// the way it is implemented now, does not always work
			if (this.modulesDir == null) {
				URL url = this.getClass().getResource("/");
				File lib = new File(url.toURI());
				logger.info("looking for " + new File(lib, TOOLS_DIST_BASENAME + "-x.y.z" + ".tar.gz"));
				File toolsTar = Files.getLatestVersion(lib, "chipster-tools", "tar.gz");
				loadFromTar(toolsTar);
			}
		
		}
	
		if (this.modulesDir != null) {
			logger.info("loading modules from " + this.modulesDir.getCanonicalPath());
		} else {
			String s = (String.format("no %s dir found after looking for: " +
					rootDir + File.separator + "%s, " + 
					rootDir + File.separator + TOOLS_DIST_BASENAME + "-x.y.z/%s, " +
					rootDir + File.separator + TOOLS_DIST_BASENAME + "-x.y.z.tar.gz" + File.separator + TOOLS_DIST_BASENAME + "-x.y.z/%s, " + 
					new File(this.getClass().getResource("/").toURI()) + File.separator + TOOLS_DIST_BASENAME + "-x.y.z" + ".tar.gz/%s",
					MODULES_DIR_NAME, MODULES_DIR_NAME, MODULES_DIR_NAME, MODULES_DIR_NAME, MODULES_DIR_NAME));
			logger.error(s);
			throw new FileNotFoundException(s);
		}
	}

	private void loadFromTar(File toolsTar) {
		if (toolsTar != null && toolsTar.exists()) {
			logger.info("found " + toolsTar);
			logger.info("looking for " + toolsTar.getAbsolutePath() + File.separator + TOOLS_DIST_BASENAME + "-x.y.z" + File.separator + MODULES_DIR_NAME);
			TFile distFile = new TFile(toolsTar);
			TFile toolsDirInTar = new TFile(Files.getLatestVersion(distFile, TOOLS_DIST_BASENAME, null));
			if (toolsDirInTar != null && toolsDirInTar.exists()) {
				File possibleModulesDir = new TFile(toolsDirInTar, MODULES_DIR_NAME);
				if (possibleModulesDir.exists() && possibleModulesDir.isDirectory()) {
					this.modulesDir = possibleModulesDir;
					logger.info("modules dir " + this.modulesDir + " found");
				} else {
					logger.info("modules dir " + possibleModulesDir + " not found");
				}
			} else {
				logger.info("no " + toolsTar.getAbsolutePath() + ":" + TOOLS_DIST_BASENAME + "-x.y.z found");
			}
		}
	}

	/**
	 * Load all the tool modules in this toolbox. Put them to the modules list.
	 * 
	 * @throws IOException
	 */
	private void loadModuleDescriptions() throws IOException {
	
		// Iterate over all module directories, and over all module files inside them
		List<String> moduleLoadSummaries = new LinkedList<String>();
		for (String moduleDirName : modulesDir.list()) {
			TFile moduleDir = new TFile(modulesDir, moduleDirName);
	
			if (moduleDir.isDirectory()) {
	
				// Load module specification files, if they exist (one module dir can contain multiple module specification files)
				for (String moduleFilename : moduleDir.list()) {
					if (moduleFilename.endsWith("-module.xml")) {
						TFile moduleFile = new TFile(moduleDir, moduleFilename);
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