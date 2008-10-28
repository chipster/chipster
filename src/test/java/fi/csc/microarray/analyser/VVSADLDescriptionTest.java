package fi.csc.microarray.analyser;

import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.LinkedList;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.analyser.r.VVSADLTool;
import fi.csc.microarray.module.chipster.ChipsterVVSADLParser.Validator;

public class VVSADLDescriptionTest {

	@BeforeSuite
	protected void setUp() throws Exception {
		MicroarrayConfiguration.loadConfiguration();
	}

	@Test(groups = {"smoke"} )
	public void testDescriptions() throws FileNotFoundException, MicroarrayException {
		
		LinkedList<String> files = new LinkedList<String>();
		files.addAll(Arrays.asList(MicroarrayConfiguration.getValues("analyser", "operations")));
		files.addAll(Arrays.asList(MicroarrayConfiguration.getValues("analyser", "hidden-operations")));
		
		for (String file : files) {
			try {
				String vvsadl;
				System.out.println(file);
				if (file.split("\\.").length > 2) {
					// class
					JavaAnalysisJobBase jobBase = (JavaAnalysisJobBase)Class.forName(file).newInstance();
					vvsadl = jobBase.getVVSADL();
				} else { 
					// script file
					if (!file.contains("/old") && !file.contains("/hidden")) {
						file = file.replace("/R", "/R/R-2.6.1"); // choose correct R script directory
					}
					VVSADLTool.ParsedRScript res = new VVSADLTool().parseRScript(getClass().getResourceAsStream(file));
					vvsadl = res.VVSADL;
				}
				new Validator().validate(file, vvsadl);
			} catch (Exception e) {
				e.printStackTrace();
				Assert.fail("when parsing " + file + ": " + e.getMessage() + " (" + e.getClass().getSimpleName() + ")");
			}
		}

	}
}
