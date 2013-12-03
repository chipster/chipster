package fi.csc.microarray.description;

import java.io.IOException;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.exception.MicroarrayException;

public class SADLTokeniserTest {

	@Test
	public void testTokenising() throws MicroarrayException, IOException {
		SADLTokeniser tokens = new SADLTokeniser("TOOL chainsaw:\"This is a \\\"chainsaw\\\" (electric)\" (Chainsaw can be used for \"cutting\" trees (like pinetrees\\).)", "");
		Assert.assertEquals(tokens.next(), "TOOL");
		Assert.assertEquals(tokens.next(), "chainsaw");
		Assert.assertEquals(tokens.next(), ":");
		Assert.assertEquals(tokens.next(), "This is a \"chainsaw\" (electric)");
		Assert.assertEquals(tokens.next(), "Chainsaw can be used for \"cutting\" trees (like pinetrees).");
	}
	

	public static void main(String[] args) throws MicroarrayException, IOException, IllegalConfigurationException {
		DirectoryLayout.initialiseSimpleLayout().getConfiguration();
		new SADLTokeniserTest().testTokenising();
	}
}
