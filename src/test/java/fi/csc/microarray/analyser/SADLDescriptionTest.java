package fi.csc.microarray.analyser;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.LinkedList;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.ChipsterSADLParser.Validator;

public class SADLDescriptionTest {

	@BeforeSuite
	protected void setUp() throws Exception {
		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"comp"}));			
	}

	@Test(groups = {"smoke"} )
	public void testDescriptions() throws FileNotFoundException, MicroarrayException {
		
		LinkedList<String> files = new LinkedList<String>();
		files.addAll(Arrays.asList(DirectoryLayout.getInstance().getConfiguration().getStrings("comp", "operations")));
		files.addAll(Arrays.asList(DirectoryLayout.getInstance().getConfiguration().getStrings("comp", "hidden-operations")));
		
		for (String file : files) {
			try {
				String sadl;
				System.out.println("validating " + file);
				if (file.split("\\.").length > 2) {
					// class
					JavaAnalysisJobBase jobBase = (JavaAnalysisJobBase)Class.forName(file).newInstance();
					sadl = jobBase.getSADL();
				} else { 
					// script file
					if (!file.contains("/old") && !file.contains("/hidden")) {
						file = file.replace("/R", "/R-2.6.1"); // choose correct R script directory
					} else {
						file = file.replace("/R", "");
					}
					SADLTool.ParsedRScript res = new SADLTool().parseRScript(getClass().getResourceAsStream(file));
					sadl = res.SADL;
				}

				new Validator().validate(file, sadl);
			} catch (Exception e) {
				e.printStackTrace();
				Assert.fail("when parsing " + file + ": " + e.getMessage() + " (" + e.getClass().getSimpleName() + ")");
			}
		}

	}
}
