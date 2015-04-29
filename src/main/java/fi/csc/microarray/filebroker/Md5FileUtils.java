package fi.csc.microarray.filebroker;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;

public class Md5FileUtils {
	
	private static final String DELIMITER = "  ";
	private static final String NEW_LINE = "\n";
	private static final int MD5_LENGTH = 32;

	/**
	 * <p>Write md5 file atomically</p>
	 * 
	 * <p>Writes a similar file what "md5sum FILENAME > FILENAME.md5" creates on command line.
	 * Atomicity is ensured by writing to a temporary file and then renaming it.
	 * After successful write, md5 file will be in the same directory with the original file and named according
	 * to original file appended by ".md5".</p> 
	 * 
	 * <p>Exception is not thrown when an identical files exists already</p> 
	 * 
	 * @param md5 checksum of the original dataFile
	 * @param dataFile original file
	 * @throws IOException writing failed
	 */
	static public void writeMd5(String md5, File dataFile) throws IOException {
		
		File md5File = getMd5File(dataFile);
		File tempFile;
		
		tempFile = File.createTempFile(md5File.getName(), null, md5File.getParentFile());

		try (BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(tempFile), "utf-8"))) {
			String line = md5 + DELIMITER + dataFile.getName() + NEW_LINE;
			writer.write(line);
			if (!tempFile.renameTo(md5File)) {
				// rename failed
				tempFile.delete();	
				// is there an identical file already?
				if (!line.equals(readMd5(dataFile))) {
					throw new IOException("rename failed");
				} 
			}
		} catch (ChecksumParseException e) {
			throw new IOException("rename and read failed", e);
		}							
	}
	
	static public String readMd5WithoutParseException(File dataFile) throws IOException {
		String md5 = null;
		try {
			md5 = readMd5(dataFile);
		} catch (ChecksumParseException e) {
		}
		
		return md5;
	}

	/**
	 * Read a .md5 file of the given data file, in a format written by writeMd5() method.
	 * 
	 * @param dataFile
	 * @return null if md5 file doesn't exist
	 * @throws FileBrokerException if the file format is not what it should be
	 * @throws IOException  
	 */
	static public String readMd5(File dataFile) throws ChecksumParseException, IOException {
		
		File md5File = getMd5File(dataFile);
		
		if (!md5File.exists()) {
			return null;
		}
						
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(md5File), "utf-8"))) {
			String line = reader.readLine(); //this eats new line character
			String[] parts = line.split(DELIMITER);
			
			if (parts.length != 2) {
				throw new ChecksumParseException("md5 file " + dataFile + " is corrupted, not two strings: " + line);				
			}
			if (parts[0].length() != MD5_LENGTH){
				throw new ChecksumParseException("md5 file " + dataFile + " is corrupted, md5 length: " + parts[0].length() + ", " + MD5_LENGTH + " expected (" + line + ")");				
			}
			
			return parts[0];							
		}
	}
	
	static File getMd5File(File file) {
		return new File(file.getAbsolutePath() + ".md5");
	}

	public static void removeMd5(File dataFile) {
		File md5File = getMd5File(dataFile);
		md5File.delete();
	}

	public static String calculateMd5(File file) throws IOException {
		
		BufferedInputStream fileStream = new BufferedInputStream(new FileInputStream(file));
			
		try (ChecksumInputStream md5Stream = new ChecksumInputStream(fileStream, true)) {
			while (md5Stream.read(new byte[1024*1024]) != -1) {				
			}
			return md5Stream.getChecksum();
		}
	}

	/**
	 * Compare two checksums and throw ChecksumException if both checksums aren't equal.
	 * Null values are always accepted.
	 * 
	 * @param checksum1
	 * @param checksum2
	 */
	public static void verify(
		String checksum1,
		String checksum2,
		String checksum3,
		Long long1,
		Long long2,
		Long long3
		) throws ChecksumException, ContentLengthException
	{
		if( checksum1 != null && checksum1.equals("MAGIC_DO_NOT_CHECK_THIS_FILE_MD5") ) {
			return;
		}
		if( checksum2 != null && checksum2.equals("MAGIC_DO_NOT_CHECK_THIS_FILE_MD5") ) {
			return;
		}
		if( checksum3 != null && checksum3.equals("MAGIC_DO_NOT_CHECK_THIS_FILE_MD5") ) {
			return;
		}
		verify(checksum1,checksum2,checksum3);
		verify(long1,long2,long3);
		return;
	}
	
	/**
	 * Compare two checksums and throw ChecksumException if both checksums aren't equal.
	 * Null values are always accepted.
	 * 
	 * @param checksum1
	 * @param checksum2
	 * @throws ChecksumException
	 */
	public static void verify(String checksum1, String checksum2) throws ChecksumException {
		if (equalsOrNull(checksum1, checksum2)) {
			return; //ok
		} else {
			throw new ChecksumException();
		}
	}
	
	/**
	 * Compare three checksums and throw ChecksumException if checksums are not equal.
	 * Null values are always accepted.
	 * 
	 * @param checksum1
	 * @param checksum2
	 * @param checksum3
	 * @throws ChecksumException
	 */
	public static void verify(String checksum1, String checksum2, String checksum3) throws ChecksumException {
		if (equalsOrNull(checksum1, checksum2, checksum3)) {
			return; //ok
		} else {
			throw new ChecksumException();
		}
	}
	
	/**
	 * Compare three checksums and return true if checksum are equal.
	 * Null values are always accepted.
	 * 
	 * @param checksum1
	 * @param checksum2
	 * @param checksum3
	 * @return
	 */
	public static boolean equalsOrNull(String checksum1, String checksum2, String checksum3) {
		return equalsOrNull(checksum1, checksum2) && equalsOrNull(checksum2, checksum3) && equalsOrNull(checksum1, checksum3);
	}
	
	/**
	 * Compare two checksums and return true if both checksums are equal.
	 * Null values are always accepted.
	 * 
	 * @param checksum1
	 * @param checksum2
	 * @return
	 */
	public static boolean equalsOrNull(String checksum1, String checksum2) {
		if (checksum1 != null && checksum2 != null) {
			return checksum1.equals(checksum2);
		}
		//can't compare
		return true;
	}
		
	/**
	 * Compare two file sizes and throw ContentLengthException if both sizes aren't equal.
	 * Null values are always accepted.
	 * 
	 * @param long1
	 * @param long2
	 * @throws ContentLengthException
	 */
	public static void verify(Long long1, Long long2) throws ContentLengthException {
		if (equalsOrNull(long1, long2)) {
			return; //ok
		} else {
			throw new ContentLengthException();
		}
	}
	
	/**
	 * Compare three file sizes and throw ContentLengthException if sizes are not equal.
	 * Null values are always accepted.
	 * 
	 * @param long1
	 * @param long2
	 * @param long3
	 * @throws ContentLengthException
	 */
	public static void verify(Long long1, Long long2, Long long3) throws ContentLengthException {
		if (equalsOrNull(long1, long2, long3)) {
			return; //ok
		} else {
			throw new ContentLengthException();
		}
	}
	
	/**
	 * Compare three file sizes and return true if sizes are equal.
	 * Null values are always accepted.
	 * 
	 * @param long1
	 * @param long2
	 * @param long3
	 * @return
	 */
	public static boolean equalsOrNull(Long long1, Long long2, Long long3) {
		return equalsOrNull(long1, long2) && equalsOrNull(long2, long3) && equalsOrNull(long1, long3);
	}
	
	/**
	 * Compare two file sizes and return true if both sizes are equal.
	 * Null values are always accepted.
	 * 
	 * @param long1
	 * @param long2
	 * @return
	 */
	public static boolean equalsOrNull(Long long1, Long long2) {
		if (long1 != null && long2 != null) {
			return long1.equals(long2);
		}
		//can't compare
		return true;
	}
}
