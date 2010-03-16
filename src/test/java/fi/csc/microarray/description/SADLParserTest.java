package fi.csc.microarray.description;

import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

public class SADLParserTest {

	@BeforeTest
	public void init() throws IOException, IllegalConfigurationException {
		DirectoryLayout.initialiseClientLayout().getConfiguration();		
	}

	@Test(groups = {"unit"} )
	public void testVVSADLCompatibility() throws MicroarrayException, IOException {
		String vvsadl = "ANALYSIS \"Test utilities\"/\"Test tool\" (Just a test analysis for development. These descriptions are sometimes very\n" + 
				"long and might get hard to read.)\n" + 
				"INPUT CDNA microarray[...].txt OUTPUT results.txt, messages.txt\n" + 
				"PARAMETER value1 INTEGER FROM 0 TO 200 DEFAULT 10 (the first value of the result set)\n" + 
				"PARAMETER value2 DECIMAL FROM 0 TO 200 DEFAULT 20 (the second value of the result set)\n" + 
				"PARAMETER value3 DECIMAL FROM 0 TO 200 DEFAULT 30.2 (the third value of the result set)\n" + 
				"PARAMETER method PERCENT DEFAULT 34 (how much we need)\n" + 
				"PARAMETER method [linear, logarithmic, exponential] DEFAULT logarithmic (which scale to use)\n" + 
				"PARAMETER genename STRING DEFAULT at_something (which gene we are interested in)\n" + 
				"PARAMETER key COLNAME (which column we use as a key)"; 

		SADLDescription parsedDescription = new ChipsterSADLParser().parse(vvsadl);
		Assert.assertNotNull(parsedDescription);
		
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
		
		SADLDescription parsedDescription = new ChipsterSADLParser().parse(sadl);
		Assert.assertNotNull(parsedDescription);
		Assert.assertEquals(parsedDescription.getName().getID(), "util-test.R");
		Assert.assertEquals(parsedDescription.getName().getDisplayName(), "Test tool");
		Assert.assertTrue(parsedDescription.getComment().startsWith("Just a test analysis"));
		Assert.assertEquals(parsedDescription.inputs().size(), 2);
		Assert.assertTrue(parsedDescription.inputs().get(0).getName().isNameSet());
		Assert.assertEquals(parsedDescription.inputs().get(1).getName().getID(), "phenodata.tsv");
		Assert.assertEquals(parsedDescription.inputs().get(1).getName().getDisplayName(), "Experiment description");
		Assert.assertEquals(parsedDescription.inputs().get(1).getType().getName(), "GENERIC");
		Assert.assertEquals(parsedDescription.outputs().size(), 2);
		Assert.assertEquals(parsedDescription.parameters().size(), 7);
		Assert.assertEquals(parsedDescription.parameters().get(4).getSelectionOptions().length, 3);
	}

	@Test(groups = {"unit"} )
	public void testRoundtrip() throws MicroarrayException, IOException {

		// create description
		SADLDescription description = new SADLDescription(Name.createName("name", "longname"), "main comment"); 
		description.addInput(new Input(ChipsterInputTypes.GENE_EXPRS, Name.createName("input1", "input1")));
		description.addInput(new Input(ChipsterInputTypes.GENE_EXPRS, Name.createNameSet("input2", ".ext", "input set 2")));
		description.addMetaInput(new Input(ChipsterInputTypes.GENE_EXPRS, Name.createName("metainput1", "metainput1")));
		description.addMetaInput(new Input(ChipsterInputTypes.GENE_EXPRS, Name.createNameSet("metainput2", ".ext", "meta input set 2")));
		description.addOutput(new Output(Name.createName("output1","output1")));
		description.addMetaOutput(new Output(Name.createName("metaoutput1", "metaoutput1")));
		description.addParameter(new Parameter(Name.createName("parameter1", "parameter1"), ParameterType.DECIMAL, null, "1", "3", "2", "param comment 1"));
		description.addParameter(new Parameter(Name.createName("parameter2", "parameter2"), ParameterType.ENUM, new String[] {"1", "2", "2"}, null, null, "2", "param comment 2"));
		
		// do some checks to created description
		Assert.assertEquals(description.inputs().get(0).getName().getID(), "input1");
		Assert.assertEquals(description.inputs().get(1).getName().getPrefix(), "input2");
		Assert.assertEquals(description.inputs().size(), 2);
		Assert.assertEquals(description.parameters().size(), 2);
		
		// serialise
		String string = description.toString();
		
		// deserialise
		SADLDescription parsedDescription = new ChipsterSADLParser().parse(string);
				
		// serialise again
		String anotherString = parsedDescription.toString();

		// compare the two serialised versions
		Assert.assertEquals(string, anotherString);
	}
	
	public static void main(String[] args) throws MicroarrayException, IOException, IllegalConfigurationException {
		DirectoryLayout.initialiseClientLayout().getConfiguration();
		new SADLParserTest().testParsing();
		new SADLParserTest().testRoundtrip();
	}
}
