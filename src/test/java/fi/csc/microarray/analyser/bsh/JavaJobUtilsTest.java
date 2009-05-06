package fi.csc.microarray.analyser.bsh;

import java.io.File;
import java.io.IOException;

import org.testng.annotations.Test;

public class JavaJobUtilsTest {

	@Test
	public void testGetGeneNames() throws IOException {

		for (String name: JavaJobUtils.getGeneNames(new File("test.tsv"))) {
			System.out.println(name);
		}
	}
}
