package fi.csc.microarray.util;

import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.Test;

public class StringsTest {
		
	@Test(groups = { "unit"})
	public void test() {
		Assert.assertEquals("004", Strings.toString(4, 3));
		Assert.assertEquals("-004", Strings.toString(-4, 3));
		Assert.assertEquals("4", Strings.toString(4, 1));
	}

	@Test(groups = { "unit"})
	public void testRemoveEmptyLinesFromBeginning() throws IOException {
		String expect = "jes on kivaa\njoo\n\n\n";
		String input = "\n\n\t\n     \n" + expect;

		Assert.assertEquals(Strings.removeEmptyLinesFromBeginning(input), expect);
	}

	


}
