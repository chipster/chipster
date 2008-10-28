package fi.csc.microarray.util;

import java.io.File;
import java.io.IOException;
import java.util.Timer;

import org.testng.annotations.Test;

public class FileCleanUpTest {

	@Test(groups = {"unit"} )
	public void testCleanUp() throws IOException {
		Timer t = new Timer(true);
		t.schedule(new FileCleanUpTimerTask(File.createTempFile("test", ""), 60*1000), 0, 1000);
		
	}


	public static void main(String[] args) throws IOException {
		Timer t = new Timer(false);
		t.schedule(new FileCleanUpTimerTask(File.createTempFile("test", ""), 60*1000), 0, 1000);

	}
}
