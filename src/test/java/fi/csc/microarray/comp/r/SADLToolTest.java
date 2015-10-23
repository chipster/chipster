package fi.csc.microarray.comp.r;

import java.io.ByteArrayInputStream;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.chipster.toolbox.SADLTool;
import fi.csc.chipster.toolbox.SADLTool.ParsedScript;

public class SADLToolTest {

	@Test
	public void testEvilContent() throws Exception {	
		// no empty line after header
		String source = "# ANALYSIS Test/Test ()\ncontent\n";
		ParsedScript parseRScript = new SADLTool("#").parseScript(new ByteArrayInputStream(source.getBytes()));
		Assert.assertEquals(parseRScript.SADL.length(), 22);
		Assert.assertEquals(parseRScript.source.length(), source.length());
	}
	
	public static void main(String[] args) throws Exception {
		new SADLToolTest().testEvilContent();
	}
}
