
package fi.csc.chipster.toolbox;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.xml.sax.SAXException;

import fi.csc.chipster.toolbox.SADLTool.ParsedScript;
import fi.csc.chipster.toolbox.toolpartsparser.ToolPartsParser;
import fi.csc.chipster.util.StringUtils;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLGenerator;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.description.SADLSyntax.ParameterType;
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
	private List<ToolboxTool> toolsWithCustomRuntime = new LinkedList<ToolboxTool>();

	private Path moduleDir;
	private Path moduleFile;

	private String summary = null;
	private int totalCount;
	private int successfullyLoadedCount;
	private int hiddenCount;
	private int disabledCount;

	private String moduleName = null;

	private SADLReplacements sadlReplacements;

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

	public ToolboxModule(Path moduleDir, Path moduleFile, File toolsBin)
			throws ParserConfigurationException, FileNotFoundException, SAXException, IOException {
		this.moduleFile = moduleFile;
		this.moduleDir = moduleDir;
		this.sadlReplacements = new SADLReplacements(toolsBin);
		load();
	}

	public ModuleDescriptionMessage getModuleDescriptionMessage() {

		// Construct description message using the current state
		ModuleDescriptionMessage msg = new ModuleDescriptionMessage(moduleName);

		for (ToolboxCategory toolboxCategory : categories) {
			Category category = new Category(toolboxCategory.getName(), toolboxCategory.getColor(),
					toolboxCategory.isHidden());

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

	public int getTotalCount() {
		return totalCount;
	}

	public int getSuccessfullyLoadedCount() {
		return successfullyLoadedCount;
	}

	public int getHiddenCount() {
		return hiddenCount;
	}

	public int getDisabledCount() {
		return disabledCount;
	}

	/**
	 * Parses a module file and loads all tools listed in it.
	 * 
	 * @param moduleDir  directory where module specification file and runtime
	 *                   specific tool directories are located in
	 * @param moduleFile module specification file
	 * @return summary of tool counts
	 * 
	 * @throws FileNotFoundException
	 * @throws SAXException
	 * @throws IOException
	 * @throws ParserConfigurationException
	 */
	private void load() throws FileNotFoundException, SAXException, IOException, ParserConfigurationException {
		// Initialise stats
		this.totalCount = 0;
		this.successfullyLoadedCount = 0;
		this.hiddenCount = 0;
		this.disabledCount = 0;

		// Load module description file
		Document document = XmlUtil.parseReader(Files.newBufferedReader(moduleFile));
		Element moduleElement = (Element) document.getElementsByTagName("module").item(0);

		// Load and check module name
		this.moduleName = moduleElement.getAttribute("name");
		if (moduleName.isEmpty()) {
			this.summary = "not loading a module without a name";
			logger.warn(summary);
			return;
		}

		// Load categories and tools in them
		for (Element categoryElement : XmlUtil.getChildElements(moduleElement, "category")) {

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
			for (Element toolElement : XmlUtil.getChildElements(categoryElement, "tool")) {
				totalCount++;

				String toolId = toolElement.getTextContent().trim();
				if (toolId.isEmpty()) {
					continue;
				}

				// Tool visibility
				boolean toolHidden = categoryHidden; // currently tools can be hidden only if their category is hidden

				// Tool module
				Path toolModuleDir;
				String nonDefaultModuleName = toolElement.getAttribute("module");
				if (nonDefaultModuleName != null && !nonDefaultModuleName.equals("")) {
					toolModuleDir = moduleDir.getParent().resolve(nonDefaultModuleName);
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
					logger.warn("not loading " + toolId + ": parameters not ok");
					continue;
				}

				// get parser

				ToolPartsParser partsParser = Toolbox.getToolPartsParser(toolId);
				if (partsParser == null) {
					logger.warn("not loading " + toolId + ": no parser");
					continue;
				}

				// parse script parts (sadl / source)
				ParsedScript parsedScript;
				try {
					parsedScript = partsParser.parse(toolModuleDir, toolId);
				} catch (Exception e) {
					logger.warn("loading " + toolId + " failed, parsing parts failed", e);
					continue;
				}

				// check SADL by parsing it, also get tool id
				SADLDescription sadlDescription;
				try {
					sadlDescription = new ChipsterSADLParser().parse(parsedScript.SADL);
				} catch (ParseException e) {
					logger.warn("loading " + toolId + " failed, parsing sadl failed", e);
					continue;
				}

				// check that sadl tool id matches module.xml tool id
				String toolIdFromSadl = sadlDescription.getName().getID();
				if (!toolId.equals(toolIdFromSadl)) {
					logger.warn("not loading " + toolId + " mismatching tool id from module file and sadl: " + toolId
							+ " " + toolIdFromSadl);
					continue;
				}

				// get runtime from sadl or use default
				String runtimeName = sadlDescription.getRuntime();
				if (runtimeName == null) {
					runtimeName = Toolbox.getDefaultRuntime(toolId);
				}
				if (runtimeName == null) {
					logger.warn("not loading " + toolId + " could not get runtime name");
				}

				// Check that filenames are unique. Overwriting input files is a bad idea when
				// the input file is
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
					logger.warn("not loading " + toolId + ": non-unique filename(s)");
					continue;
				}

				// apply replacements (i.e. generate genome lists)
				try {
					for (Parameter param : sadlDescription.getParameters()) {
						if (param.getType() == ParameterType.ENUM) {

							List<Name> options = Arrays.asList(param.getSelectionOptions());
							param.setSelectionOptions(sadlReplacements.processNames(options));

							List<String> defaults = Arrays.asList(param.getDefaultValues());
							param.setDefaultValues(sadlReplacements.processStrings(defaults));
						}
					}
				} catch (IOException e) {
					logger.warn("not loading " + toolId + ": " + e.getMessage());
					continue;
				}

				// check parameter defaults, needs to happen after replacements
				if (!checkParameterDefaults(sadlDescription)) {
					logger.warn("not loading " + toolId + ": invalid parameter default(s)");
					continue;
				}

				String generatedSadl = null;
				try {
					// generate the SADL with the replacements (original available in
					// parsedScript.SADL)
					generatedSadl = SADLGenerator.generate(sadlDescription);
				} catch (Exception e) {
					logger.warn("not loading " + toolId + ": ", e);
					continue;
				}

				// Register the tool, override existing
				ToolboxTool toolboxTool = new ToolboxTool(toolId, sadlDescription, generatedSadl, parsedScript.code,
						parsedScript.source, moduleDir.getFileName().toString(), runtimeName);
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

				if (sadlDescription.getRuntime() != null) {
					toolsWithCustomRuntime.add(toolboxTool);
				}

				logger.debug(String.format("loaded %s %s from %s %s", toolId,
						sadlDescription.getName().getDisplayName(), toolId, hiddenStatus));
			}

		}

		// Update summary
		this.summary = "loaded " + moduleName + " " + successfullyLoadedCount + "/" + totalCount + " tools, "
				+ disabledCount + " disabled, " + hiddenCount + " hidden";
		logger.info(summary);
		if (!toolsWithCustomRuntime.isEmpty()) {
			for (ToolboxTool tool : toolsWithCustomRuntime) {
				logger.info("custom runtime for " + tool.getId() + " [" + tool.getRuntime() + "]");
			}
		}

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

	private boolean checkParameterDefaults(SADLDescription sadl) {
		boolean valid = true;
		for (Parameter param : sadl.getParameters()) {
			// parameters with default value(s)
			if (param.getDefaultValues() != null && param.getDefaultValues().length > 0) {
				switch (param.getType()) {
				case ENUM:
					// get valid options
					Set<String> options = Arrays.stream(param.getSelectionOptions()).map(name -> name.getID())
							.collect(Collectors.toSet());

					// check if default values are included in the valid options
					for (String defaultValue : param.getDefaultValues()) {
						if (!options.contains(defaultValue)) {
							valid = false;
							logger.warn(sadl.getName().getID() + " parameter " + param.getName().getID()
									+ " has invalid default " + defaultValue + " ["
									+ options.stream().collect(Collectors.joining(", ")) + "]");
						}
					}
					break;
				case INTEGER:
					break;
				default:
					break;
				}
			}
		}

		return valid;
	}

}