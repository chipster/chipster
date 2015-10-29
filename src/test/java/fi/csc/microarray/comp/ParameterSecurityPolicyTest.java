package fi.csc.microarray.comp;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.comp.r.RCompJob;

public class ParameterSecurityPolicyTest {

	@Test
	public void test() throws Exception {
		Assert.assertFalse("stop()".matches(RCompJob.RParameterSecurityPolicy.NUMERIC_VALUE_PATTERN));
		Assert.assertFalse("\"stop()".matches(RCompJob.RParameterSecurityPolicy.TEXT_VALUE_PATTERN));
		Assert.assertTrue("hgu133ahsentrezg(hgu133a)".matches(RCompJob.RParameterSecurityPolicy.TEXT_VALUE_PATTERN));

		// BeanShell and Java analysis jobs have no patterns to test
	}
	
	public static void main(String[] args) throws Exception {
		new ParameterSecurityPolicyTest().test();
		System.out.println("ParameterSecurityPolicyTest OK");
	}
}
