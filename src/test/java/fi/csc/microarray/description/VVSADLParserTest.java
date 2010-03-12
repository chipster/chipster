package fi.csc.microarray.description;

import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;
import fi.csc.microarray.module.chipster.ChipsterVVSADLParser;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.AnnotatedName;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

public class VVSADLParserTest {

	@BeforeTest
	public void init() throws IOException, IllegalConfigurationException {
		DirectoryLayout.initialiseClientLayout().getConfiguration();		
	}

	@Test(groups = {"unit"} )
	public void testParsing() throws MicroarrayException, IOException {
		String sadl = "TOOL util-test.R: \"Test tool\" (Just a test analysis for development. These descriptions are sometimes very\n" + 
				"long and might get hard to read.)\n" + 
				"INPUT microarray{...}.tsv: \"Raw data files\" TYPE CDNA\n" + 
				"INPUT phenodata.tsv: \"Experiment description\" TYPE GENERIC\n" + 
				"OUTPUT result{...}.txt\n" + 
				"OUTPUT OPTIONAL error.txt\n" + 
				"PARAMETER value1: \"The first value\" TYPE INTEGER FROM 0 TO 200 DEFAULT 10 (the first value of the result set)\n" + 
				"PARAMETER value2: \"The second value\" TYPE DECIMAL FROM 0 TO 200 DEFAULT 20 (the second value of the result set)\n" + 
				"PARAMETER OPTIONAL value3: \"The third value\" TYPE DECIMAL FROM 0 TO 200 DEFAULT 30.2 (the third value of the result set)\n" + 
				"PARAMETER method: \"The fourth value\" TYPE PERCENT DEFAULT 34 (how much we need)\n" + 
				"PARAMETER method: \"The method\" TYPE [linear, logarithmic: \"Logarithmic scale\", exponential: \"Exponential scale\"] FROM 1 TO 2 DEFAULT logarithmic (which scale to use)\n" + 
				"PARAMETER genename: \"Gene name\" TYPE STRING DEFAULT at_something (which gene we are interested in)\n" + 
				"PARAMETER key: \"Key column\" TYPE COLUMN_SEL (which column we use as a key)"; 
		
		SADLDescription parsedDescription = new ChipsterVVSADLParser().parse(sadl);
		Assert.assertNotNull(parsedDescription);
		Assert.assertEquals(parsedDescription.getAnnotatedName().getName(), "util-test.R");
		Assert.assertEquals(parsedDescription.getAnnotatedName().getHumanReadableName(), "Test tool");
		Assert.assertTrue(parsedDescription.getComment().startsWith("Just a test analysis"));
		Assert.assertEquals(parsedDescription.inputs().size(), 2);
		Assert.assertTrue(parsedDescription.inputs().get(0).getAnnotatedName().isNameSet());
		Assert.assertEquals(parsedDescription.inputs().get(1).getAnnotatedName().getName(), "phenodata.tsv");
		Assert.assertEquals(parsedDescription.inputs().get(1).getAnnotatedName().getHumanReadableName(), "Experiment description");
		Assert.assertEquals(parsedDescription.inputs().get(1).getType().getName(), "GENERIC");
		Assert.assertEquals(parsedDescription.outputs().size(), 2);
		Assert.assertEquals(parsedDescription.parameters().size(), 7);
		Assert.assertEquals(parsedDescription.parameters().get(4).getSelectionOptions().length, 3);
	}

	@Test(groups = {"unit"} )
	public void testRoundtrip() throws MicroarrayException, IOException {

		// create description
		SADLDescription description = new SADLDescription(AnnotatedName.createName("name", "longname"), "main comment"); 
		description.addInput(new Input(ChipsterInputTypes.GENE_EXPRS, AnnotatedName.createName("input1", "input1")));
		description.addInput(new Input(ChipsterInputTypes.GENE_EXPRS, AnnotatedName.createNameSet("input2", ".ext", "input set 2")));
		description.addMetaInput(new Input(ChipsterInputTypes.GENE_EXPRS, AnnotatedName.createName("metainput1", "metainput1")));
		description.addMetaInput(new Input(ChipsterInputTypes.GENE_EXPRS, AnnotatedName.createNameSet("metainput2", ".ext", "meta input set 2")));
		description.addOutput("output1");
		description.addMetaOutput("metaoutput1");
		description.addParameter(new Parameter(AnnotatedName.createName("parameter1", "parameter1"), ParameterType.DECIMAL, null, "1", "3", "2", "param comment 1"));
		description.addParameter(new Parameter(AnnotatedName.createName("parameter2", "parameter2"), ParameterType.ENUM, new String[] {"1", "2", "2"}, null, null, "2", "param comment 2"));
		
		// do some checks to created description
		Assert.assertEquals(description.inputs().get(0).getAnnotatedName().getName(), "input1");
		Assert.assertEquals(description.inputs().get(1).getAnnotatedName().getPrefix(), "input2");
		Assert.assertEquals(description.inputs().size(), 2);
		Assert.assertEquals(description.parameters().size(), 2);
		
		// serialise
		String string = description.toString();
		
		// deserialise
		SADLDescription parsedDescription = new ChipsterVVSADLParser().parse(string);
				
		// serialise again
		String anotherString = parsedDescription.toString();

		// compare the two serialised versions
		Assert.assertEquals(string, anotherString);
	}
	
	public static void main(String[] args) throws MicroarrayException, IOException, IllegalConfigurationException {
		DirectoryLayout.initialiseClientLayout().getConfiguration();
		new VVSADLParserTest().testParsing();
		new VVSADLParserTest().testRoundtrip();
	}
}
