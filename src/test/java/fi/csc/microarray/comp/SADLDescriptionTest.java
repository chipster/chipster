package fi.csc.microarray.comp;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.junit.Assert;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import fi.csc.chipster.toolbox.SADLTool;
import fi.csc.chipster.toolbox.SADLTool.ParsedScript;
import fi.csc.microarray.comp.java.JavaCompJobBase;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.module.chipster.ChipsterSADLParser.Validator;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.XmlUtil;

public class SADLDescriptionTest {

	public static void main(String[] args) throws Exception {
		SADLDescriptionTest test = new SADLDescriptionTest();
		test.testDescriptions();
		test.testSADLTool();
		System.out.println("SADLDescriptionTest OK");
	}

	@Test
	public void testSADLTool() throws IOException {
		
		String rScript = 
			"# TOOL \"BLAST\" / blastn.sadl: blastn (Heuristic tool to search hits for a nucleotide sequence from a nucleotode sequence database.)" + "\n" +
			"# INPUT query: \"Query sequence\" TYPE GENERIC" + "\n" +
			"# OUTPUT out.txt" + "\n" +
			"# PARAMETER OPTIONAL query_loc: \"Location on the query sequence\" TYPE STRING (Location of the search region on the query sequence. Format: start-stop, for example: 23-66.  Default: the whole query sequence)" + "\n" +
			"# PARAMETER OPTIONAL strand: \"Query strand\" TYPE [both: Both, minus: Minus, plus: Plus] DEFAULT both ( Query strand or strands to search against the database.  Default: both strands.)" + "\n" +
			"" + "\n" +
			"echo('foo');" + "\n";
		
		SADLTool tool = new SADLTool("#");
		ParsedScript parsedScript = tool.parseScript(new ByteArrayInputStream(rScript.getBytes()));
		String rScript2 = tool.toScriptString(parsedScript);
		
		Assert.assertEquals(rScript2, rScript);

	}

	@Test
	public void testDescriptions() throws Exception {

		DirectoryLayout.uninitialise();
		DirectoryLayout.initialiseUnitTestLayout();
		
		class ToolSpec {
			private String module;
			private String resource;
			private String runtime;
			private String runtimeDir;
			private String toolSpecificModule;
			private boolean isHidden;
			
			public ToolSpec(String module, String resource, String runtime, String runtimeDir, 
					String toolSpecificModule, boolean isHidden) {
				this.module = module;
				this.resource = resource;
				this.runtime = runtime;
				this.runtimeDir = runtimeDir;
				this.toolSpecificModule = toolSpecificModule;
				this.isHidden = isHidden;
			}
		}
		
		LinkedList<ToolSpec> toolSpecs = new LinkedList<ToolSpec>();
		

		// Find all files to check
		// We don't use ToolRepository to do all this because it assumes different directory layout

		// Iterate over runtimes to get their dirs
		Map<String, String> runtimeDirMap = new HashMap<String, String>();
		Document runtimeDoc = XmlUtil.parseFile(new File("src/main/applications/wrapper/comp/conf/runtimes.xml"));
		NodeList runtimes = runtimeDoc.getDocumentElement().getElementsByTagName("runtime");
		for (int i = 0; i < runtimes.getLength(); i++) {
			Element runtime = (Element)runtimes.item(i);
			String name = runtime.getElementsByTagName("name").item(0).getTextContent();
			String dir = null;
			NodeList parameters = runtime.getElementsByTagName("parameter");
			for (int j = 0; j < parameters.getLength(); j++) {
				if ("toolPath".equals(((Element)parameters.item(j)).getElementsByTagName("name").item(0).getFirstChild().getTextContent())) {
					dir = ((Element)parameters.item(j)).getElementsByTagName("value").item(0).getFirstChild().getTextContent();
					break;
				}
			}
			runtimeDirMap.put(name, dir);		
		}
		
		// Iterate through all tools and collect their resource definitions (filename, classname etc.) 
		for (File file : Files.listFilesRecursively(new File("src/main/modules/"))) {
			if (file.getName().endsWith("-module.xml")) {
				Document module = XmlUtil.parseFile(file);
				String moduleName = module.getDocumentElement().getAttribute("name");
				NodeList categories = module.getDocumentElement().getElementsByTagName("category");
				for (int j = 0; j < categories.getLength(); j++) {
					Element category = (Element)categories.item(j);
					boolean isHidden = "true".equals(category.getAttribute("hidden"));
					NodeList tools = category.getElementsByTagName("tool");
					for (int i = 0; i < tools.getLength(); i++) {
						Element tool = (Element)tools.item(i);
						String runtimeName = tool.getAttribute("runtime");
						String moduleOverride = tool.getAttribute("module");
						String toolSpecificModule = moduleOverride.isEmpty() ? moduleName : moduleOverride;
						String resource = tool.getElementsByTagName("resource").item(0).getTextContent();
						toolSpecs.add(new ToolSpec(moduleName, resource.trim(), runtimeName, runtimeDirMap.get(runtimeName), toolSpecificModule, isHidden));
					}
				}
			}
		}
		
		// Load SADL from each resource
		String missingManuals = "";
		for (ToolSpec toolspec : toolSpecs) {
			try {
				String sadl = null;
				
				// Try to figure out what it was. We have this information in the module descriptions,
				// but for keeping the test code simple, we infer it from the resource name. Currently
				// that is enough, in future we might need to use the actual module loading facility
				// to parse the module files.
				boolean isFile;
				String manualName;
				if ("java".equals(toolspec.runtime)) {
					// Is a class name

					System.out.println("validating class " + toolspec.resource + " in " + toolspec.module);
					JavaCompJobBase jobBase = (JavaCompJobBase)Class.forName(toolspec.resource).newInstance();
					sadl = jobBase.getSADL();
					isFile = false;
					manualName = toolspec.resource.substring(toolspec.resource.lastIndexOf('.') + 1);
					
				} else { 
					// Is a file name
					isFile = true;
					manualName = toolspec.resource.substring(0, toolspec.resource.indexOf('.'));
					
					// Determine which file it is
					File file = new File(new File("src/main/modules"), toolspec.toolSpecificModule + File.separator + toolspec.runtimeDir + File.separator + toolspec.resource);

					// Determine file type and process it
					if (file.getName().endsWith(".R") || file.getName().endsWith(".py")) {
						// Is an R script
						SADLTool.ParsedScript res = new SADLTool("#").parseScript(new FileInputStream(file));

						sadl = res.SADL;

					} else if (file.getName().endsWith(".bsh")) {
						// Is a BeanShell script
						SADLTool.ParsedScript res = new SADLTool("//").parseScript(new FileInputStream(file));
						sadl = res.SADL;

					} else if (file.getName().endsWith(".acd")) {
						continue; // nothing to validate in machine generated descriptions (EMBOSS)

					} else {
						// IS a plain SADL file
						sadl = Files.fileToString(file);
					}

					System.out.println("validating file " + file.getCanonicalFile() + " in " + toolspec.module);
				}

				// Finally, validate the description
				if (sadl != null) {
					List<SADLDescription> descriptions = new Validator().validate(toolspec.resource, sadl);
					Assert.assertEquals(1, descriptions.size());
					if (isFile) {
						Assert.assertEquals(toolspec.resource, descriptions.get(0).getName().getID());
					}
					if (!toolspec.isHidden && !new File("src/main/manual/" + manualName + ".html").exists()) {
						missingManuals += "\nManual page missing for " + toolspec.resource;
					}
					
				} else {
					throw new RuntimeException("don't know what to do with: " + toolspec);
				}
				
			} catch (Exception e) {
				e.printStackTrace();
				Assert.fail("when parsing " + toolspec + ": " + e.getMessage() + " (" + e.getClass().getSimpleName() + ")");
			}
		}
		System.out.println(missingManuals);

	}
}
