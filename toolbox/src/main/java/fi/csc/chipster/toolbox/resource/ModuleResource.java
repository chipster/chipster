package fi.csc.chipster.toolbox.resource;

import java.io.IOException;
import java.util.LinkedHashSet;

import javax.inject.Singleton;
import javax.ws.rs.GET;
import javax.ws.rs.Path;
import javax.ws.rs.Produces;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;

import com.fasterxml.jackson.databind.JsonNode;
import com.fasterxml.jackson.databind.node.ArrayNode;
import com.fasterxml.jackson.databind.node.JsonNodeFactory;
import com.fasterxml.jackson.databind.node.ObjectNode;

import fi.csc.chipster.toolbox.Toolbox;
import fi.csc.chipster.toolbox.ToolboxModule;
import fi.csc.chipster.toolbox.ToolboxModule.ToolboxCategory;
import fi.csc.chipster.toolbox.rest.RestUtils;
import fi.csc.chipster.toolbox.ToolboxService;
import fi.csc.chipster.toolbox.ToolboxTool;

@Singleton
@Path("modules")
public class ModuleResource {

	private Toolbox toolbox;

	public ModuleResource(Toolbox toolbox) throws IOException {
		this.toolbox = toolbox;
	}

	@GET
//    @JacksonFeatures(serializationEnable =  { SerializationFeature.INDENT_OUTPUT })
	@Produces(MediaType.APPLICATION_JSON)
	public final Response getModules() throws IOException {

		// used several times in this method
		// use local reference to avoid toolbox update messing things up
		Toolbox localToolbox = this.toolbox;

		// node factory for creating json nodes
		JsonNodeFactory factory = new JsonNodeFactory(false);

		ArrayNode modules = factory.arrayNode();

		// modules

		// define order or the modules
		String[] priorityModules = new String[] { "ngs", "microarray", "misc" };
		LinkedHashSet<ToolboxModule> orderedModules = new LinkedHashSet<ToolboxModule>();
		for (String name : priorityModules) {
			ToolboxModule module = localToolbox.getModule(name);
			if (module != null) {
				orderedModules.add(module);
			}
		}
		for (ToolboxModule module : localToolbox.getModules()) {
			if (!orderedModules.contains(module)) {
				orderedModules.add(module);
			}
		}

		// go through modules
		for (ToolboxModule toolboxModule : orderedModules) {
			ObjectNode module = factory.objectNode();
			module.put("name", toolboxModule.getNamePretty());

			// categories in a module
			ArrayNode categories = factory.arrayNode();
			module.set("categories", categories);
			for (ToolboxCategory toolboxCategory : toolboxModule.getCategories()) {
				ObjectNode category = factory.objectNode();
				category.put("name", toolboxCategory.getName());
				category.put("color", toolboxCategory.getColor());
				category.put("hidden", toolboxCategory.isHidden());

				// tools in a category
				ArrayNode tools = factory.arrayNode();
				category.set("tools", tools);
				for (ToolboxTool toolboxTool : toolboxCategory.getTools()) {
					JsonNode tool = RestUtils.getObjectMapper().convertValue(toolboxTool.getSadlDescription(),
							JsonNode.class);
					tools.add(tool);
				}
				categories.add(category);
			}

			modules.add(module);
		}

		return Response.ok(modules).build();
	}

	@GET
	@Path("zip")
	@Produces(MediaType.APPLICATION_OCTET_STREAM)
	public Response getZip() {
		return Response.ok(toolbox.getZipStream())
				// hint filename
				.header("Content-Disposition", "attachment; filename=\"" + ToolboxService.TOOLS_ZIP_NAME + "\"")
				.build();

	}

	public void setToolbox(Toolbox newToolbox) {
		this.toolbox = newToolbox;
	}

}
