package fi.csc.microarray.description;

import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLSyntax.ParameterType;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class SADLParserTest {

	@Test
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
				"PARAMETER key COLUMN_SEL (which column we use as a key)"; 

		SADLDescription parsedDescription = new ChipsterSADLParser().parse(vvsadl);
		Assert.assertNotNull(parsedDescription);
		Assert.assertEquals(parsedDescription.getName().getID(), "Test tool");
		Assert.assertEquals(parsedDescription.getName().getDisplayName(), "Test tool");
		Assert.assertTrue(parsedDescription.getDescription().startsWith("Just a test analysis"));
		Assert.assertEquals(parsedDescription.getInputs().size(), 1);
		Assert.assertTrue(parsedDescription.getInputs().get(0).getName().isNameSet());
		Assert.assertEquals(parsedDescription.getOutputs().size(), 2);
		Assert.assertEquals(parsedDescription.getParameters().size(), 7);
		Assert.assertEquals(parsedDescription.getParameters().get(4).getSelectionOptions().length, 3);
		
	}

	@Test
	public void testParsing() throws MicroarrayException, IOException {
		String sadl = "TOOL util-test.R: \"Test tool\" (Just a test analysis for development. These descriptions are sometimes very\n" + 
				"long and might get hard to read. (Note that certain operators must be escaped.\\))\n" + 
				"INPUT microarray{...}.tsv: \"Raw data files\" TYPE CDNA\n" + 
				"INPUT phenodata.tsv: \"Experiment description\" TYPE GENERIC\n" + 
				"OUTPUT result{...}.txt: \"Result files\"\n" + 
				"OUTPUT OPTIONAL error.txt: \"Error, if any\" (Here's an example of a comment on an output file.)\n" + 
				"PARAMETER value1: \"The first value\" TYPE INTEGER FROM 0 TO 200 DEFAULT 10 (the first value of the result set)\n" + 
				"PARAMETER value2: \"The second value\" TYPE DECIMAL FROM 0 TO 200 DEFAULT 20 (the second value of the result set)\n" + 
				"PARAMETER OPTIONAL value3: \"The third value\" TYPE DECIMAL FROM 0 TO 200 DEFAULT 30.2 (the third value of the result set)\n" + 
				"PARAMETER value4: \"The fourth value\" TYPE PERCENT DEFAULT 34 (how much we need)\n" + 
				"PARAMETER method: \"The enumeration\" TYPE [option1: \"First option\", option2: \"Second option\", option3: \"Third option\"] FROM 1 TO 2 DEFAULT option1, option2 (which options are selected)\n" + 
				"PARAMETER genename: \"Gene name\" TYPE STRING DEFAULT at_something (which gene we are interested in)\n" + 
				"PARAMETER key: \"Key column\" TYPE COLUMN_SEL (which column we use as a key)"; 
		
		// do a comprehensive test of the parsed description
		SADLDescription parsedDescription = new ChipsterSADLParser().parse(sadl);

		// TOOL
		Assert.assertNotNull(parsedDescription);
		Assert.assertEquals(parsedDescription.getName().getID(), "util-test.R"); // test name parts here, not repeated after this
		Assert.assertEquals(parsedDescription.getName().getDisplayName(), "Test tool"); // test name parts here, not repeated after this
		Assert.assertEquals(parsedDescription.getName().toString(), "util-test.R: \"Test tool\"");
		Assert.assertTrue(parsedDescription.getDescription().startsWith("Just a test analysis"));
		Assert.assertTrue(parsedDescription.getDescription().endsWith("must be escaped.)"));
		
		// INPUTS
		Assert.assertEquals(parsedDescription.getInputs().size(), 2);
		Assert.assertTrue(parsedDescription.getInputs().get(0).getName().isNameSet());
		Assert.assertEquals(parsedDescription.getInputs().get(0).getName().getPrefix(), "microarray");
		Assert.assertEquals(parsedDescription.getInputs().get(0).getName().getPostfix(), ".tsv");
		Assert.assertEquals(parsedDescription.getInputs().get(0).getName().getDisplayName(), "Raw data files");
		Assert.assertEquals(parsedDescription.getInputs().get(0).getType().getName(), "CDNA");
		Assert.assertFalse(parsedDescription.getInputs().get(0).isOptional());
		Assert.assertEquals(parsedDescription.getInputs().get(1).getName().toString(), "phenodata.tsv: \"Experiment description\"");
		Assert.assertEquals(parsedDescription.getInputs().get(1).getType().getName(), "GENERIC");
		Assert.assertFalse(parsedDescription.getInputs().get(1).isOptional());
		
		// OUTPUTS
		Assert.assertEquals(parsedDescription.getOutputs().size(), 2);
		Assert.assertTrue(parsedDescription.getOutputs().get(0).getName().isNameSet());
		Assert.assertEquals(parsedDescription.getOutputs().get(0).getName().getPrefix(), "result");
		Assert.assertEquals(parsedDescription.getOutputs().get(0).getName().getPostfix(), ".txt");
		Assert.assertEquals(parsedDescription.getOutputs().get(0).getName().getDisplayName(), "Result files");
		Assert.assertFalse(parsedDescription.getOutputs().get(0).isOptional());
		Assert.assertEquals(parsedDescription.getOutputs().get(1).getName().toString(), "error.txt: \"Error, if any\"");
		Assert.assertEquals(parsedDescription.getOutputs().get(1).getDescription(), "Here's an example of a comment on an output file.");
		Assert.assertTrue(parsedDescription.getOutputs().get(1).isOptional());
		
		// PARAMETERS
		Assert.assertEquals(parsedDescription.getParameters().size(), 7);
		Assert.assertEquals(parsedDescription.getParameters().get(0).getName().toString(), "value1: \"The first value\"");
		Assert.assertEquals(parsedDescription.getParameters().get(0).getType(), ParameterType.INTEGER);
		Assert.assertEquals(parsedDescription.getParameters().get(0).getFrom(), "0");
		Assert.assertEquals(parsedDescription.getParameters().get(0).getTo(), "200");
		Assert.assertEquals(parsedDescription.getParameters().get(0).getDefaultValue(), "10");
		Assert.assertEquals(parsedDescription.getParameters().get(0).getDescription(), "the first value of the result set");
		Assert.assertFalse(parsedDescription.getParameters().get(0).isOptional());
		Assert.assertEquals(parsedDescription.getParameters().get(1).getName().toString(), "value2: \"The second value\"");
		Assert.assertEquals(parsedDescription.getParameters().get(1).getType(), ParameterType.DECIMAL);
		Assert.assertEquals(parsedDescription.getParameters().get(1).getFrom(), "0");
		Assert.assertEquals(parsedDescription.getParameters().get(1).getTo(), "200");
		Assert.assertEquals(parsedDescription.getParameters().get(1).getDefaultValue(), "20");
		Assert.assertEquals(parsedDescription.getParameters().get(1).getDescription(), "the second value of the result set");
		Assert.assertFalse(parsedDescription.getParameters().get(1).isOptional());
		Assert.assertEquals(parsedDescription.getParameters().get(2).getName().toString(), "value3: \"The third value\"");
		Assert.assertEquals(parsedDescription.getParameters().get(2).getType(), ParameterType.DECIMAL);
		Assert.assertEquals(parsedDescription.getParameters().get(2).getFrom(), "0");
		Assert.assertEquals(parsedDescription.getParameters().get(2).getTo(), "200");
		Assert.assertEquals(parsedDescription.getParameters().get(2).getDefaultValue(), "30.2");
		Assert.assertEquals(parsedDescription.getParameters().get(2).getDescription(), "the third value of the result set");
		Assert.assertTrue(parsedDescription.getParameters().get(2).isOptional());
		Assert.assertEquals(parsedDescription.getParameters().get(3).getName().toString(), "value4: \"The fourth value\"");
		Assert.assertEquals(parsedDescription.getParameters().get(3).getType(), ParameterType.PERCENT);
		Assert.assertEquals(parsedDescription.getParameters().get(3).getDefaultValue(), "34");
		Assert.assertEquals(parsedDescription.getParameters().get(3).getDescription(), "how much we need");
		Assert.assertFalse(parsedDescription.getParameters().get(3).isOptional());
		Assert.assertEquals(parsedDescription.getParameters().get(4).getName().toString(), "method: \"The enumeration\"");
		Assert.assertEquals(parsedDescription.getParameters().get(4).getType(), ParameterType.ENUM);
		Assert.assertEquals(parsedDescription.getParameters().get(4).getDefaultValues().length, 2);
		Assert.assertEquals(parsedDescription.getParameters().get(4).getDefaultValues()[0], "option1");
		Assert.assertEquals(parsedDescription.getParameters().get(4).getDefaultValues()[1], "option2");
		Assert.assertEquals(parsedDescription.getParameters().get(4).getDescription(), "which options are selected");
		Assert.assertEquals(parsedDescription.getParameters().get(4).getFrom(), "1");
		Assert.assertEquals(parsedDescription.getParameters().get(4).getTo(), "2");
		Assert.assertEquals(parsedDescription.getParameters().get(4).getSelectionOptions().length, 3);
		Assert.assertEquals(parsedDescription.getParameters().get(4).getSelectionOptions()[0].toString(), "option1: \"First option\"");
		Assert.assertEquals(parsedDescription.getParameters().get(4).getSelectionOptions()[1].toString(), "option2: \"Second option\"");
		Assert.assertEquals(parsedDescription.getParameters().get(4).getSelectionOptions()[2].toString(), "option3: \"Third option\"");
		Assert.assertFalse(parsedDescription.getParameters().get(4).isOptional());
		Assert.assertEquals(parsedDescription.getParameters().get(5).getName().toString(), "genename: \"Gene name\"");
		Assert.assertEquals(parsedDescription.getParameters().get(5).getType(), ParameterType.STRING);
		Assert.assertEquals(parsedDescription.getParameters().get(5).getDefaultValue(), "at_something");
		Assert.assertEquals(parsedDescription.getParameters().get(5).getDescription(), "which gene we are interested in");
		Assert.assertFalse(parsedDescription.getParameters().get(5).isOptional());
		Assert.assertEquals(parsedDescription.getParameters().get(6).getName().toString(), "key: \"Key column\"");
		Assert.assertEquals(parsedDescription.getParameters().get(6).getType(), ParameterType.COLUMN_SEL);
		Assert.assertEquals(parsedDescription.getParameters().get(6).getDescription(), "which column we use as a key");
		Assert.assertFalse(parsedDescription.getParameters().get(6).isOptional());
	}

	@Test
	public void testRoundtrip() throws MicroarrayException, IOException {

		// create description
		SADLDescription description = generateDescription();
		
		// serialise
		String string = description.toString();

		// deserialise
		SADLDescription parsedDescription = new ChipsterSADLParser().parse(string);
				
		// serialise again
		String anotherString = parsedDescription.toString();

		// compare the two serialised versions
		Assert.assertEquals(string, anotherString);
		
		// to guard against serialisation omissions, we should do complete check between description and parsedDescription
		// now we just do some checks
		Assert.assertEquals(parsedDescription.getInputs().get(0).getName().getID(), "input1");
		Assert.assertEquals(parsedDescription.getInputs().get(0).getDescription(), "input comment");
		Assert.assertTrue(parsedDescription.getInputs().get(0).isOptional());
		Assert.assertEquals(parsedDescription.getInputs().get(1).getName().getPrefix(), "input2");
		Assert.assertEquals(parsedDescription.getInputs().size(), 4);
		Assert.assertEquals(parsedDescription.getParameters().size(), 3);
		Assert.assertEquals(parsedDescription.getParameters().get(1).getFrom(), "1");
		Assert.assertEquals(parsedDescription.getParameters().get(2).getDefaultValue(), "");
		
	}

	@Test
	public void testEscapes() throws MicroarrayException, IOException {

		// create description and check
		SADLDescription description = generateDescription();		
		Assert.assertEquals(description.getDescription(), "main comment (funny)");
		
		// serialise and check
		String string = description.toString();
		Assert.assertTrue(string.contains("main comment (funny\\)"));
		
		// deserialise and check
		SADLDescription parsedDescription = new ChipsterSADLParser().parse(string);
		Assert.assertEquals(parsedDescription.getDescription(), "main comment (funny)");
				
		// serialise again and check
		String anotherString = parsedDescription.toString();
		Assert.assertTrue(anotherString.contains("main comment (funny\\)"));
	}
	
	private SADLDescription generateDescription() {
		SADLDescription description = new SADLDescription(Name.createName("name", "longname/displayname"), "main comment (funny)");
		Input input = new Input(ChipsterInputTypes.GENE_EXPRS, Name.createName("input1", "input1"), true);
		input.setDescription("input comment");
		description.addInput(input);
		description.addInput(new Input(ChipsterInputTypes.GENE_EXPRS, Name.createNameSet("input2", ".ext", "input set 2")));
		description.addInput(new Input(ChipsterInputTypes.GENE_EXPRS, Name.createName("metainput1", "metainput1"), false, true));
		description.addInput(new Input(ChipsterInputTypes.GENE_EXPRS, Name.createNameSet("metainput2", ".ext", "meta input set 2"), false, true));
		description.addOutput(new Output(Name.createName("output1","output1")));
		description.addOutput(new Output(Name.createName("metaoutput1", "metaoutput1"), false, true));
		description.addParameter(new Parameter(Name.createName("parameter1", "parameter1"), ParameterType.DECIMAL, null, "1", "3", "2", "param comment 1"));
		description.addParameter(new Parameter(Name.createName("parameter2", "parameter2"), ParameterType.ENUM, new Name[] {Name.createName("1"), Name.createName("2"), Name.createName("3")}, "1", "2", new String[]{"1", "2"}, "param comment 2"));
		description.addParameter(new Parameter(Name.createName("parameter3", "parameter3"), ParameterType.STRING, null, null, null, "", "empty default value"));
		return description;
	}
	
	public static void main(String[] args) throws MicroarrayException, IOException, IllegalConfigurationException {
		
		new SADLParserTest().testParsing();
		new SADLParserTest().testRoundtrip();
		new SADLParserTest().testEscapes();
		new SADLParserTest().testVVSADLCompatibility();
		System.out.println("SADLParserTest OK");
	}
}
