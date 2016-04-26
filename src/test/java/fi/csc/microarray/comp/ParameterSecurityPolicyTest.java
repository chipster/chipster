package fi.csc.microarray.comp;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.comp.python.PythonCompJob;
import fi.csc.microarray.comp.r.RCompJob;
/**
 * BeanShell and Java analysis jobs have no patterns to test
 * 
 */
public class ParameterSecurityPolicyTest {

	@Test
	public void testText() throws Exception {

		String [] reject = { 
				"\"stop()",
				"\"",
				"'",
				"\t",
				"\n",
				"$",
				"{}",
				"\u2600",
				};
		
		String [] accept = { 
				"abc",
				"123",
				"hgu133ahsentrezg(hgu133a)",
				"jou()",
				"*",
				"+",
				",",
				".",
				"-",
				"_",
				" ",
				"åÖä",
				"é",
				"",
				"\u0073",
				"\u00E3",
				};

		
		
		for (String s : reject) {
			System.out.print("testing " + s);
			Assert.assertFalse(s.matches(RCompJob.RParameterSecurityPolicy.TEXT_VALUE_PATTERN));
			Assert.assertFalse(s.matches(PythonCompJob.PythonParameterSecurityPolicy.TEXT_VALUE_PATTERN));
			System.out.println (" -> rejected");
		}

		for (String s : accept) {
			System.out.print("testing " + s);
			Assert.assertTrue(s.matches(RCompJob.RParameterSecurityPolicy.TEXT_VALUE_PATTERN));
			Assert.assertTrue(s.matches(PythonCompJob.PythonParameterSecurityPolicy.TEXT_VALUE_PATTERN));
			System.out.println (" -> accepted");
		}
	}


	@Test
	public void testNumber() throws Exception {

		String [] reject = { 
				"\"stop()",
				"\"",
				"'",
				"abc",
				"$",
				"{}",
				"\u2600",
				};
		
		String [] accept = { 
				"0",
				"123",
				"-0",
				"-10000000000000000000000000000000",
				"1.7",
				"1324234.7009808",
				};

		
		
		for (String s : reject) {
			System.out.print("testing " + s);
			Assert.assertFalse(s.matches(RCompJob.RParameterSecurityPolicy.NUMERIC_VALUE_PATTERN));
			System.out.println (" -> rejected");
		}

		for (String s : accept) {
			System.out.print("testing " + s);
			Assert.assertTrue(s.matches(RCompJob.RParameterSecurityPolicy.NUMERIC_VALUE_PATTERN));
			System.out.println (" -> accepted");
		}
	}

	
	
	public static void main(String[] args) throws Exception {
		ParameterSecurityPolicyTest t = new ParameterSecurityPolicyTest(); 
		t.testText();
		t.testNumber();
		System.out.println("ParameterSecurityPolicyTest OK");
	}
}
