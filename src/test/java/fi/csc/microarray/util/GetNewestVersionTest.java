package fi.csc.microarray.util;

import java.io.File;

import org.junit.Test;

public class GetNewestVersionTest {

	@Test
	public void testGetNewestVersion() {
		File f = Files.getNewestVersion(new File("/tmp/neppi"), "chipster-tools", "tar");
		System.out.println(f.getName());
	}

	@Test
	public void testVersion() {
		System.out.println(Version.parse("1.2.33"));
	
		System.out.println(Version.parse("2.22.1").compareTo(Version.parse("2.1")));
	
	}
}
