package fi.csc.microarray.client.dataimport;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.client.dataimport.trimmer.ConditionalStringReplace;
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
	
	@Test(groups = {"unit"} )
	public void testConditional() {
		DataTrimmer trimmer = new DataTrimmer();
		
		trimmer.pushOperation(new ConditionalStringReplace(Double.MIN_VALUE, 0.4, false, false, "P", 1));
		trimmer.pushOperation(new ConditionalStringReplace(0.4, 0.8, true, false, "M", 1));
		trimmer.pushOperation(new ConditionalStringReplace(0.8, Double.MAX_VALUE, true, false, "A", 1));
		
		String column1 = "0.1";
		String column2 = "0.3";
		String column3 = "0.4";
		String column4 = "0.5";
		String column5 = "0.6";
		String column6 = "0.8";
		String column7 = "0.9";
		
		String trimmedColumn1 = trimmer.doTrimming(column1, 1);
		String trimmedColumn2 = trimmer.doTrimming(column2, 1);
		String trimmedColumn3 = trimmer.doTrimming(column3, 1);
		String trimmedColumn4 = trimmer.doTrimming(column4, 1);
		String trimmedColumn5 = trimmer.doTrimming(column5, 1);
		String trimmedColumn6 = trimmer.doTrimming(column6, 1);
		String trimmedColumn7 = trimmer.doTrimming(column7, 1);
		
		System.out.println("Conditional replacement");
		System.out.println("Column1 ("+column1+"): " + trimmedColumn1);
		System.out.println("Column2 ("+column2+"): " + trimmedColumn2);
		System.out.println("Column3 ("+column3+"): " + trimmedColumn3);
		System.out.println("Column4 ("+column4+"): " + trimmedColumn4);
		System.out.println("Column5 ("+column5+"): " + trimmedColumn5);
		System.out.println("Column6 ("+column6+"): " + trimmedColumn6);
		System.out.println("Column7 ("+column7+"): " + trimmedColumn7);
		
		Assert.assertTrue(trimmedColumn1.equals("P"));
		Assert.assertTrue(trimmedColumn2.equals("P"));
		Assert.assertTrue(trimmedColumn3.equals("M"));
		Assert.assertTrue(trimmedColumn4.equals("M"));
		Assert.assertTrue(trimmedColumn5.equals("M"));
		Assert.assertTrue(trimmedColumn6.equals("A"));
		Assert.assertTrue(trimmedColumn7.equals("A"));
	}
	
}
