package fi.csc.microarray.analyser;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.analyser.emboss.EmbossAnalysisJob;
import fi.csc.microarray.analyser.r.RAnalysisJob;
import fi.csc.microarray.analyser.shell.ShellAnalysisJob;

public class ParameterSecurityPolicyTest {

	@Test
	public void test() throws Exception {
		Assert.assertFalse("stop()".matches(RAnalysisJob.RParameterSecurityPolicy.NUMERIC_VALUE_PATTERN));
		Assert.assertFalse("\"stop()".matches(RAnalysisJob.RParameterSecurityPolicy.TEXT_VALUE_PATTERN));
		Assert.assertTrue("hgu133ahsentrezg(hgu133a)".matches(RAnalysisJob.RParameterSecurityPolicy.TEXT_VALUE_PATTERN));
		Assert.assertFalse(";ls".matches(ShellAnalysisJob.ShellParameterSecurityPolicy.COMMAND_LINE_SAFE_VALUE_PATTERN));
		Assert.assertFalse(";ls".matches(EmbossAnalysisJob.EmbossParameterSecurityPolicy.COMMAND_LINE_SAFE_VALUE_PATTERN));
		Assert.assertTrue("embl:X65923".matches(EmbossAnalysisJob.EmbossParameterSecurityPolicy.COMMAND_LINE_SAFE_VALUE_PATTERN));
		Assert.assertTrue("OS=Homo sapiens".matches(EmbossAnalysisJob.EmbossParameterSecurityPolicy.COMMAND_LINE_SAFE_VALUE_PATTERN));

		// BeanShell and Java analysis jobs have no patterns to test
	}
	
	public static void main(String[] args) throws Exception {
		new ParameterSecurityPolicyTest().test();
		System.out.println("ParameterSecurityPolicyTest OK");
	}
}
