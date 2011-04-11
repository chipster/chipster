package fi.csc.microarray.analyser;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.LinkedList;

import org.testng.Assert;
import org.testng.annotations.Test;
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

	@Test(groups = {"unit"} )
	private void testSADLTool() throws IOException {
		
		String rScript = 
			"# TOOL \"BLAST\" / blastn.sadl: blastn (Heuristic tool to search hits for a nucleotide sequence from a nucleotode sequence database.)" + "\n" +
			"# INPUT query: \"Query sequence\" TYPE GENERIC" + "\n" +
			"# OUTPUT out.txt" + "\n" +
			"# PARAMETER OPTIONAL query_loc: \"Location on the query sequence\" TYPE STRING (Location of the search region on the query sequence. Format: start-stop, for example: 23-66.  Default: the whole query sequence)" + "\n" +
			"# PARAMETER OPTIONAL strand: \"Query strand\" TYPE [both: Both, minus: Minus, plus: Plus] DEFAULT both ( Query strand or strands to search against the database.  Default: both strands.)" + "\n" +
			"" + "\n" +
			"echo('jee');" + "\n";
		
		SADLTool tool = new SADLTool("#");
		ParsedScript parsedScript = tool.parseScript(new ByteArrayInputStream(rScript.getBytes()));
		String rScript2 = tool.toScriptString(parsedScript);
		
		Assert.assertEquals(rScript2, rScript);

	}

	@Test(groups = {"unit"} )
	public void testDescriptions() throws Exception {
				
		DirectoryLayout.initialiseUnitTestLayout();
		
		LinkedList<String> resources = new LinkedList<String>();
		
		// Iterate through all tools and collect their resource definitions (filename, classname etc.) 
		for (File file : new File("src/main/applications/wrapper/comp/conf").listFiles()) {
			if (file.getName().endsWith("-module.xml")) {
				Document module = XmlUtil.parseFile(file);
				NodeList tools = module.getDocumentElement().getElementsByTagName("tool");
				for (int i = 0; i < tools.getLength(); i++) {
					Element tool = (Element)tools.item(i);
					String resource = tool.getElementsByTagName("resource").item(0).getTextContent();
					
					if (resource.endsWith(".acd")) {
						// Refers to EMBOSS ACD, they are converted in the code level, so there is nothing to test
						continue;
					}

					resources.add(resource);
				}
			}
		}
		
		// Load SADL from each resource
		for (String resource : resources) {
			try {
				String sadl = null;
				
				// Try to figure out what it was. We have this information in the module descriptions,
				// but for keeping the test code simple, we infer it from the resource name. Currently
				// that is enough, in future we might need to use the actual module loading facility
				// to parse the module files.
				if (resource.split("\\.").length > 2) {
					// Is a class name

					System.out.println("validating class " + resource);
					JavaAnalysisJobBase jobBase = (JavaAnalysisJobBase)Class.forName(resource).newInstance();
					sadl = jobBase.getSADL();
					
				} else { 
					// Is a file name
					
					// Collect all possible files that the resource name might refer to
					File[] dirsContainingDescriptions = new File[] {
						new File("src/main/modules/chipster/bsh"),
						new File("src/main/modules/sequence/shell"),
						new File("src/main/modules/chipster/R-2.12"),
						new File("src/main/modules/ngs/R-2.12"),
						new File("src/main/modules/chipster/R-2.10"),
					};
					LinkedList<File> potentialFiles = new LinkedList<File>(); 

					for (File dir : dirsContainingDescriptions) {
						for (File file : dir.listFiles()) {
							potentialFiles.add(file);
						}
					}
					
					// Determine which file it is
					for (File file : potentialFiles) {
						if (file.getName().endsWith(resource)) {

							// Found the file, determine the type and process it
							if (resource.endsWith(".R")) {
								// Is an R script
								SADLTool.ParsedScript res = new SADLTool("#").parseScript(new FileInputStream(file));
								
								sadl = res.SADL;

							} else if (resource.endsWith(".bsh")) {
								// Is a BeanShell script
								SADLTool.ParsedScript res = new SADLTool("//").parseScript(new FileInputStream(file));
								sadl = res.SADL;

							} else {
								// IS a plain SADL file
								sadl = Files.fileToString(file);
							}

							System.out.println("validating file " + file.getCanonicalFile());
							
							break; // we are done with this file
						}
					}
					
				}

				// Finally, validate the description
				if (sadl != null) {
					new Validator().validate(resource, sadl);
					
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
