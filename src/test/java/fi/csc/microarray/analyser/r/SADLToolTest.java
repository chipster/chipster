package fi.csc.microarray.analyser.r;

import java.io.ByteArrayInputStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.SADLTool;
import fi.csc.microarray.analyser.SADLTool.ParsedScript;

public class SADLToolTest {

	@Test(groups = { "unit"})
	public void testEvilContent() throws Exception {	
		// no empty line after header
		String source = "# ANALYSIS Test/Test ()\ncontent\n";
		ParsedScript parseRScript = new SADLTool().parseScript(new ByteArrayInputStream(source.getBytes()), "#");
		Assert.assertEquals(parseRScript.SADL.length(), 22);
		Assert.assertEquals(parseRScript.source.length(), source.length());
	}
	
	public static void main(String[] args) throws Exception {
		new SADLToolTest().testEvilContent();
	}
}
