package fi.csc.microarray.util;

import java.io.File;

import org.junit.Test;

import de.schlichtherle.truezip.file.TFile;

public class GetLatestVersionTest {

	@Test
	public void testGetLatestVersion() {
		File f = Files.getLatestVersion(new File("/tmp/neppi"), "chipster-tools", "tar.gz");
		System.out.println(f.getName());
	}

	@Test
	public void testGetLatestVersionDir() {
		File f = Files.getLatestVersion(new File("/tmp/neppi"), "chipster-tools", null);
		System.out.println(f.getName());
	}

	@Test
	public void testGetLatestVersionDirFromTar() {
		File f = Files.getLatestVersion(new TFile("/Users/hupponen/git/chipster/chipster-tools-3.6.3.tar.gz"), "chipster-tools", null);
		System.out.println(f.getName());
	}

	
	@Test
	public void testVersion() {
		System.out.println(Version.parse("1.2.33"));
		System.out.println(Version.parse("2.22.1").compareTo(Version.parse("2.1")));
	}
}
