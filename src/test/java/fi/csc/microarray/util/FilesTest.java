package fi.csc.microarray.util;

import java.io.File;
import java.io.IOException;
import java.util.Calendar;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.AfterTest;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

public class FilesTest {

	private File testRoot, dir_1, dir_2;
	
	@BeforeTest
	public void setUp() throws IOException {
		// yes, yes, not perfect
		testRoot = File.createTempFile("deltree-test-", "");
		testRoot.delete();
		Assert.assertTrue(testRoot.mkdir());
		System.out.println("created test root dir: " + testRoot.getAbsolutePath());
		
		dir_1 = createDirWithFiles(testRoot, "dir_1", 2);
		dir_2 = createDirWithFiles(testRoot, "dir_2", 3);
		
		Files.createSymbolicLink(new File(dir_1, "file1"), new File(dir_1, "fileLink"));
		Files.createSymbolicLink(dir_2, new File(dir_1, "dirLink"));
	}
	
	@AfterTest
	public void tearDown() {
		testRoot.delete();
	}

	
	@Test
	public void testListFilesRecursivelySortByDateOldestFirst() throws IOException, InterruptedException {
		File dir = new File("testdir");

		try {
			dir.mkdir();

			File f1 = new File(dir, "f1");
			f1.createNewFile();
			File f2 = new File(dir, "f2");
			f2.createNewFile();
			File f3 = new File(dir, "f3");
			f3.createNewFile();
			Calendar c = Calendar.getInstance(); 
			c.add(Calendar.DATE, 1);
			f2.setLastModified(c.getTime().getTime()); // f2 is newer (created "tomorrow")
			c.add(Calendar.DATE, -2);
			f3.setLastModified(c.getTime().getTime()); // f3 is older (created "yesterday")

			List<File> list = Files.listFilesRecursivelySortByDateOldestFirst(dir);
			Assert.assertEquals(list.get(0), f3);
			Assert.assertEquals(list.get(1), f1);
			Assert.assertEquals(list.get(2), f2);
			Assert.assertEquals(list.size(), 3);
			
		} finally {

			Files.delTree(dir);
		}
	}

	@Test
	public void testDelTree() throws IOException {
		
		// basic dir
		File dir = new File(testRoot, "dir");
		dir.mkdir();
		Assert.assertTrue(Files.delTree(dir));
		Assert.assertFalse(dir.exists());
		
		// basic file
		File file = new File(testRoot, "file");
		file.createNewFile();
		Assert.assertTrue(Files.delTree(file));
		Assert.assertFalse(file.exists());
	
		// mix of things
		Assert.assertTrue(Files.delTree(dir_1));
		Assert.assertFalse(dir_1.exists());
		Assert.assertTrue(Files.delTree(dir_2));
		Assert.assertFalse(dir_2.exists());
	}

	
	private File createDirWithFiles(File parentDir, String dirName, int depth) throws IOException {
		File dir = new File(parentDir, dirName);
		dir.mkdir();
		File file_1 = new File(dir, "file_1");
		file_1.createNewFile();
		File file_2 = new File(dir, "file_2");
		file_2.createNewFile();

		if (depth > 0) {
			createDirWithFiles(dir, "dir", --depth);
		}
		
		return dir;
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		new FilesTest().testListFilesRecursivelySortByDateOldestFirst();
//		FilesTest t = new FilesTest();
//		try {
//			t.setUp();
//			t.testDelTree();
//		} finally {
//			t.tearDown();
//		}
	}
}
