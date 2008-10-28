package fi.csc.microarray.util;

import org.testng.Assert;
import org.testng.annotations.Test;

public class StringsTest {
		
	@Test(groups = { "unit"})
	public void test() {
		Assert.assertEquals("004", Strings.toString(4, 3));
		Assert.assertEquals("-004", Strings.toString(-4, 3));
		Assert.assertEquals("4", Strings.toString(4, 1));
	}
}
