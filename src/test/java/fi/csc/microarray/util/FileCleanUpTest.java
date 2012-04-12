package fi.csc.microarray.util;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import java.util.Timer;

import org.apache.commons.io.FileUtils;
import org.testng.annotations.Test;

public class FileCleanUpTest {

	@Test(groups = {"unit"} )
	public void testCleanUp() throws IOException {
		Timer t = new Timer(true);
		t.schedule(new FileCleanUpTimerTask(File.createTempFile("test", ""), 60*1000), 0, 1000);
		
	}
	

	@Test(groups = {"unit"} )
	public void listFilesSortByDate() throws IOException {
		List<File> files = Files.listFilesRecursivelySortByDateOldestFirst(new File("/Users/taavi/Downloads"));
		for (File file : files) {
			System.out.println(file.getName() + ", " + FileUtils.byteCountToDisplaySize(file.length()) + ", " + new SimpleDateFormat("yyyy-MM-dd HH:mm:SS").format(new Date(file.lastModified())));
		}
	}

	


	public static void main(String[] args) throws IOException {
		new FileCleanUpTest().listFilesSortByDate();
		
	}
}
