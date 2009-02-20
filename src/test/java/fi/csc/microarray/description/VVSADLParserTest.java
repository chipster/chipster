package fi.csc.microarray.description;

import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.config.MicroarrayConfiguration;
import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;
import fi.csc.microarray.description.ParsedVVSADL.Input;
import fi.csc.microarray.description.ParsedVVSADL.Parameter;
import fi.csc.microarray.description.VVSADLSyntax.ParameterType;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;
import fi.csc.microarray.module.chipster.ChipsterVVSADLParser;

public class VVSADLParserTest {

	@BeforeTest
	public void init() throws IOException, OldConfigurationFormatException {
		MicroarrayConfiguration.loadConfiguration();		
	}
	
	@Test(groups = {"unit"} )
	public void testRoundtrip() throws MicroarrayException, IOException {

		// create description
		ParsedVVSADL description = new ParsedVVSADL("name", "package", "main comment"); 
		description.addInput(Input.createInput(ChipsterInputTypes.GENE_EXPRS, "input1"));
		description.addInput(Input.createInputSet(ChipsterInputTypes.GENE_EXPRS, "input2", ".ext"));
		description.addMetaInput(Input.createInput(ChipsterInputTypes.GENE_EXPRS, "metainput1"));
		description.addMetaInput(Input.createInputSet(ChipsterInputTypes.GENE_EXPRS, "metainput2", ".ext"));
		description.addOutput("output1");
		description.addMetaOutput("metaoutput1");
		description.addParameter(new Parameter("parameter1", ParameterType.DECIMAL, null, "1", "3", "2", "param comment 1"));
		description.addParameter(new Parameter("parameter2", ParameterType.ENUM, new String[] {"1", "2", "2"}, null, null, "2", "param comment 2"));
		
		// do some checks to created description
		Assert.assertEquals(description.inputs().get(0).getName(), "input1");
		Assert.assertEquals(description.inputs().get(1).getPrefix(), "input2");
		Assert.assertEquals(description.inputs().size(), 2);
		Assert.assertEquals(description.parameters().size(), 2);
		
		// serialise
		String string = description.toString();
		
		// deserialise
		ParsedVVSADL parsedDescription = new ChipsterVVSADLParser().parse(string);
		
		// serialise again
		String anotherString = parsedDescription.toString();
		System.out.println(anotherString);
		
		// compare the two serialised versions
		Assert.assertEquals(string, anotherString);
	}
	
	public static void main(String[] args) throws MicroarrayException, IOException, OldConfigurationFormatException {
		MicroarrayConfiguration.loadConfiguration();
		new VVSADLParserTest().testRoundtrip();
	}
}
