package fi.csc.microarray.client.dataimport;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.client.dataimport.trimmer.DataTrimmer;
import fi.csc.microarray.client.dataimport.trimmer.NormalStringReplace;
import fi.csc.microarray.client.dataimport.trimmer.ReqularExpressionStringReplace;

public class DataTrimmerTest {
	
	@Test(groups = {"unit"} ) 
	public void testNormal(){
		DataTrimmer trimmer = new DataTrimmer();
		
		trimmer.pushOperation(new NormalStringReplace("a", "A", false, 1));
		trimmer.pushOperation(new NormalStringReplace("(", "[", false, 1));
		trimmer.pushOperation(new NormalStringReplace(")", "}", false, 2));
		trimmer.pushOperation(new NormalStringReplace("A", "s", false, 1));
		trimmer.popOperation();
		trimmer.pushOperation(new NormalStringReplace("A", "f", false, 1));
		
		String column1 = "(aAaA)";
		String column2 = "(aAaA)";
		
		String trimmedColumn1 = trimmer.doTrimming(column1, 1);
		String trimmedColumn2 = trimmer.doTrimming(column2, 2);
		
		System.out.println("Normal string replacement ("+trimmer.getOperationCount()+" operations)");
		System.out.println("Column1: " + trimmedColumn1);
		System.out.println("Column2: " + trimmedColumn2);
		
		Assert.assertTrue(trimmedColumn1.equals("[ffff)"));
		Assert.assertTrue(trimmedColumn2.equals("(aAaA}"));
	}

	@Test(groups = {"unit"} )
	public void testRegexp(){
		DataTrimmer trimmer = new DataTrimmer();
		
		trimmer.pushOperation(new ReqularExpressionStringReplace("go*", "aa", 1));
		trimmer.pushOperation(new ReqularExpressionStringReplace("go+g", "a", 2));
		
		String column1 = "gogogooogo";
		String column2 = "goooooggogoo";
		
		String trimmedColumn1 = trimmer.doTrimming(column1, 1);
		String trimmedColumn2 = trimmer.doTrimming(column2, 2);
		
		System.out.println("Regular expression replacement");
		System.out.println("Column1: " + trimmedColumn1);
		System.out.println("Column2: " + trimmedColumn2);
		
		Assert.assertTrue(trimmedColumn1.equals("aaaaaaaa"));
		Assert.assertTrue(trimmedColumn2.equals("aaoo"));
	}
}
