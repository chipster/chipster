package fi.csc.microarray.analyser;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.junit.Assert;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;

import fi.csc.microarray.analyser.SADLTool.ParsedScript;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.config.DirectoryLayout;
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
		
		LinkedList<Pair<String, String>> resources = new LinkedList<Pair<String, String>>();

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
				NodeList tools = module.getDocumentElement().getElementsByTagName("tool");
				for (int i = 0; i < tools.getLength(); i++) {
					Element tool = (Element)tools.item(i);
					String runtimeName = tool.getAttribute("runtime");
					String moduleOverride = tool.getAttribute("module");
					String toolSpecificModule = moduleOverride.isEmpty() ? moduleName : moduleOverride;
					String resource = tool.getElementsByTagName("resource").item(0).getTextContent();
					String dir = runtimeDirMap.get(runtimeName);
					if (dir == null) {
						resources.add(new ImmutablePair<String, String>(runtimeName, resource.trim()));
					} else {
						resources.add(new ImmutablePair<String, String>(runtimeName, toolSpecificModule + File.separator + dir + File.separator + resource));
					}
				}
			}
		}
		
		// Load SADL from each resource
		for (Pair<String, String> resource : resources) {
			try {
				String sadl = null;
				
				// Try to figure out what it was. We have this information in the module descriptions,
				// but for keeping the test code simple, we infer it from the resource name. Currently
				// that is enough, in future we might need to use the actual module loading facility
				// to parse the module files.
				if ("java".equals(resource.getLeft())) {
					// Is a class name

					System.out.println("validating class " + resource.getRight());
					JavaAnalysisJobBase jobBase = (JavaAnalysisJobBase)Class.forName(resource.getRight()).newInstance();
					sadl = jobBase.getSADL();
					
				} else { 
					// Is a file name
					
					// Determine which file it is
					File file = new File(new File("src/main/modules"), resource.getRight());

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

					System.out.println("validating file " + file.getCanonicalFile());
				}

				// Finally, validate the description
				if (sadl != null) {
					new Validator().validate(resource.getRight(), sadl);
					
				} else {
					throw new RuntimeException("don't know what to do with: " + resource);
				}
				
			} catch (Exception e) {
				e.printStackTrace();
				Assert.fail("when parsing " + resource + ": " + e.getMessage() + " (" + e.getClass().getSimpleName() + ")");
			}
		}

	}
}
