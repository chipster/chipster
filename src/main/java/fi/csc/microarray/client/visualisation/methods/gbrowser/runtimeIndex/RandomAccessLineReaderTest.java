package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

/**
 * Test for {@link RandomAccessLineReader}. Generates a temp file and checks that RandomAccessLineReader
 * reads that file correctly.
 * 
 * @author klemela
 */
public class RandomAccessLineReaderTest {
	
	private static final int TEST_FILE_ROWS = 1000;
		
	public static void main (String[] args) throws IOException, URISyntaxException, GBrowserException {
		
		//Test empty file
		File emptyFile = File.createTempFile("RandomAccesLineReaderTest", ".txt");
		DataUrl dataUrl = new DataUrl(emptyFile);
		RandomAccessLineReader lineReader = new RandomAccessLineReader(dataUrl);
		System.out.println(lineReader.setPosition(0) == false);		
		emptyFile.delete();
		
		//Create some artificial content
		File testFile = getTestFile();	
		testFile(testFile, true);
		
		//Create some artificial content with really long rows
		testFile = getTestFile(true, 100);	
		testFile(testFile, false);	
	}

	private static void testFile(File testFile, boolean thoroughTest) throws FileNotFoundException,
			IOException, URISyntaxException, MalformedURLException,
			GBrowserException {
		RandomAccessLineReader lineReader;
		//Read the file with standard tools for reference
		List<String> lines = getTestReferenceList(testFile);
		DataUrl dataUrl = new DataUrl(testFile);
		lineReader = new RandomAccessLineReader(dataUrl);
		
		//Check that lineReader reports correct file length
		System.out.println(lineReader.length() == testFile.length());		
		
		//Check that incorrect positions are recognized 
		System.out.println(lineReader.setPosition(-1) == false);	
		System.out.println(lineReader.setPosition(testFile.length()) == false);
		
		boolean readLineException = false;
		try {
			//This should produce exception, because last setPostion was incorrect
			lineReader.readLine();
		} catch (IOException e) {
			readLineException = true;
		}
		
		System.out.println(readLineException);				

		//This is valid position
		System.out.println(lineReader.setPosition(0) == true);

		
		// Read through the whole file and check that every line is correct
		lineReader.setPosition(0);
		
		for (String refernceLine : lines) {
			String line = lineReader.readLine();
			testLineComparison(refernceLine, line);
		}

		if (thoroughTest) {

			// Find out a location in the middle of the file for buffer tests
			int testLinePosition = 0;

			for (int i = 0; i < 500; i++) {
				testLinePosition += lines.get(i).length() + 1;
			}

			//Right answers
			String[] line500 = {"M499-", "499-", "99-", "9-", "-", "" };

			//Read something little bit before the actual place to make sure that buffer is ready
			for (int i = 0; i < 6; i++) {
				int j = testLinePosition - 6 + i;

				//Fill buffer by reading little bit before the actual place
				lineReader.setPosition(j - 1024);
				lineReader.readLine();

				//This should give the last characters of line number 499

				//Read from buffer
				lineReader.setPosition(j);			
				String line = lineReader.readLine();			
				//System.out.println(line);			
				testLineComparison(line500[i], line);
			}

			//Force buffer update		
			for (int i = 0; i < 6; i++) {
				int j = testLinePosition - 6 + i;

				//Make buffer useless by reading somwhere else
				lineReader.setPosition(RandomAccessLineReader.HTTP_BUFFER_SIZE*2);
				lineReader.readLine();

				//Refresh buffer and read
				lineReader.setPosition(j);
				String line = lineReader.readLine();			
				//System.out.println(line);			
				testLineComparison(line500[i], line);
			}

			//Buffer runs out	
			for (int i = 0; i < 6; i++) {
				int j = testLinePosition - 6 + i;

				//Fill buffer so that it has only one byte of the actual content
				final int BUFFER_BYTES_LEFT = 1;
				//Fill buffer
				lineReader.setPosition(j - RandomAccessLineReader.HTTP_BUFFER_SIZE + BUFFER_BYTES_LEFT);			
				lineReader.readLine();			

				if (lineReader.setPosition(j)) {

					//Check the buffer refill works
					String line = lineReader.readLine();			
					//System.out.println(line);
					testLineComparison(line500[i], line);

				} else {
					System.err.println("setPosition failed");
				}
			}

			//Now some testing at the end of file
			//Correct answers for end of file tests
			String[] endOfFile = {"-M999-", "M999-", "999-", "99-", "9-", "-", "" };

			//End of file		
			for (int i =  0; i < 7; i++) {
				int j = i - 7;
				lineReader.setPosition(testFile.length() + j);
				String line = lineReader.readLine();			
				//System.out.println(line);
				testLineComparison(endOfFile[i], line);
			}

			System.out.println(lineReader.readLine() == null);
		}
		
		testFile.delete();
	}

	private static void testLineComparison(String refernceLine, String line) {
		if (!refernceLine.equals(line)) {
			System.err.println("Unequal lines, \n\tRef:  \t" + refernceLine + "\n\tTest: \t" + line);
		}
	}

	public static List<String> getTestReferenceList(File testFile)
			throws FileNotFoundException, IOException {
		List<String> lines = new LinkedList<String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(testFile));

		//Create reference list with standard tools 
		String line; 
		
		while ((line = reader.readLine()) != null) {
			lines.add(line);
		}			
		reader.close();

		return lines;
	}
	
	public static File getTestFile() throws IOException {
		return getTestFile(false, TEST_FILE_ROWS);
	}

	/**
	 * Create a temp file with 1000 rows:
	 * 
	 * A000-B000-C000-D000-E000-F000-G000-H000-I000-J000-K000-L000-M000-
	 * A001-B001-C001-D001-E001-F001-G001-H001-I001-J001-K001-L001-M001-
	 * A002-B002-C002-D002-E002-F002-G002-H002-I002-J002-K002-L002-M002-
	 * ...
	 * A997-B997-C997-D997-E997-F997-G997-H997-I997-J997-K997-L997-M997-
	 * A998-B998-C998-D998-E998-F998-G998-H998-I998-J998-K998-L998-M998-
	 * A999-B999-C999-D999-E999-F999-G999-H999-I999-J999-K999-L999-M999-
	 * @param testFileRows 
	 * 
	 * @return
	 * @throws IOException
	 */
	public static File getTestFile(boolean longRows, int testFileRows) throws IOException {
		//Generate test file
		File testFile = File.createTempFile("RandomAccessLineReader-test-file", ".txt");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(testFile));
		
		String[] cols = new String[] { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M" };
		
		int rowRepeat = 1;
		
		if (longRows) {
			rowRepeat = 100;
		}
		

		for (int i = 0; i < testFileRows; i++) {
			String row = String.format("%03d", i);

			String line = "";
			for (int j = 0; j < rowRepeat; j++) {
				for (String col : cols) {
					line += col + row + "-";				
				}
			}

			//System.out.println(line);
			writer.write(line);
			writer.newLine();
		}
		writer.flush();
		writer.close();
		return testFile;
	}
}
