package fi.csc.microarray.analyser.r;

import java.io.ByteArrayInputStream;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.VVSADLTool;
import fi.csc.microarray.analyser.VVSADLTool.ParsedRScript;

public class VVSADLToolTest {

	@Test(groups = { "unit"})
	public void testEvilContent() throws Exception {	
		// no empty line after header
		String source = "# ANALYSIS Test/Test ()\ncontent\n";
		ParsedRScript parseRScript = new VVSADLTool().parseRScript(new ByteArrayInputStream(source.getBytes()));
		Assert.assertEquals(parseRScript.VVSADL.length(), 22);
		Assert.assertEquals(parseRScript.rSource.length(), source.length());
	}
	
	public static void main(String[] args) throws Exception {
		new VVSADLToolTest().testEvilContent();
	}
}
