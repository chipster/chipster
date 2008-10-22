/*
 * Created on Feb 14, 2005
 *
 */
package fi.csc.microarray.util;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileFilter;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import org.apache.commons.io.filefilter.AgeFileFilter;
import org.apache.commons.io.filefilter.DirectoryFileFilter;
import org.apache.commons.io.filefilter.IOFileFilter;
import org.apache.commons.io.filefilter.OrFileFilter;
import org.mortbay.util.IO;

/**
 * @author Taavi Hupponen, Aleksi Kallio
 *
 */
public class Files {

	
	/**
	 * Deletes a file or a directory recursively.
	 * @param dir directory or file to be deleted
	 * @return true if deleting was successful, false file does not exist or deleting it failed
	 */
	public static boolean delTree(File dir) {
		try {
			// if dir is directory, make it empty recursively
			if (dir.isDirectory()) {
				File[] contents = dir.listFiles();
				if (contents.length > 0) {
					for (int i = 0; i < contents.length; i++) {
						delTree(contents[i]);
					}
				}
			}
			
			// now dir should be either a file or an empty directory
			// try to delete it
			if (dir.delete()) {
				return true;
			} else {
				return false;
			}
	
		} catch (SecurityException se) {
			return false;
		}
	}

	public static byte[] fileToBytes(File file) throws IOException  {
		if (file == null) {
			throw new IllegalArgumentException("parameter file is null.");
		}
		
		InputStream s = new FileInputStream(file);
		try {
			return inputStreamToBytes(new FileInputStream(file));	
		} finally {
			s.close();
		}
	}

	public static byte[] inputStreamToBytes(InputStream input) throws IOException  {
		if (input == null) {
			throw new NullPointerException("parameter input is null.");
		}
		
		BufferedInputStream in = new BufferedInputStream(input);
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		IO.copy(in, out);
		return out.toByteArray();
	}
	
	public static String fileToString(File file) throws IOException {
		if (file == null) {
			throw new IllegalArgumentException("parameter file is null.");
		}
		
		InputStream s = new FileInputStream(file);
		try {
			return inputStreamToString(new FileInputStream(file));	
		} finally {
			s.close();
		}
	}

	
	public static String inputStreamToString(InputStream input) throws IOException  {
		if (input == null) {
			throw new NullPointerException("parameter input is null.");
		}

		StringBuffer buffer = new StringBuffer();
		BufferedReader inputReader = new BufferedReader(new InputStreamReader(input));
		String line;
		for (line = inputReader.readLine(); line != null; line = inputReader.readLine()) {
			buffer.append(line + "\n");
		}
		
		return buffer.toString();
	}

	public static boolean sweepAndCreateDirectory(File dir) {
		
		if (dir.exists()) {
			if (!dir.isDirectory()) {
				return false;
			}
			
			if (!delTree(dir)) {
				return false;
			}
		}
		
		if (!dir.mkdir()) {
			return false;
		}

		return true;
		
	}
	
	public static boolean equalInputStreamContent(InputStream input1, InputStream input2) throws IOException {
		InputStream in1 = new BufferedInputStream(input1);
		InputStream in2  = new BufferedInputStream(input2);
		
		boolean equal = true;
		for (int byte1 = in1.read(); byte1 != -1 && equal == true; byte1 = in1.read() ) {
			if (byte1 != in2.read()) {
				equal = false;
			}
		}
		if (in2.read() != -1) {
			equal = false;
		}
		
		return equal;
		
	}
	
	/**
	 * Walks the baseDir recursively and deleted files and directories older than cutoff.
	 * 
	 * If a directory is old but contains files (which are not too old), it is not deleted.
	 * 
	 * TODO better problem handling?
	 * 
	 * @param baseDir
	 * @param cutoff milliseconds 
	 */
	public static void cleanOldFiles(File baseDir, long cutoff ) {
		walkAndDelete(baseDir, new AgeFileFilter(System.currentTimeMillis() - cutoff));
	}

	public static void walkAndDelete(File baseDir, IOFileFilter filter) {
		
		
		File[] files = baseDir.listFiles((FileFilter)new OrFileFilter(DirectoryFileFilter.INSTANCE, filter));
		
		if (files == null) {
			return;
		}
		
		for (File f: files) {
			if (f.isDirectory()) {
				
				// check the filter for directory before walking it as walking might affect the filter
				// conditions such as directory modified time
				boolean toBeDeleted = false;
				if (filter.accept(f)) {
					toBeDeleted = true;
				}

				// walk into dir
				walkAndDelete(f, filter);
				
				// possibly delete dir 
				if (toBeDeleted) {
					f.delete();
				}
			
			} else {
				f.delete();
			}
		}
	}

	/**
	 * Parses the whole filename and returns name and extension separately.
	 * Returns array with possibly empty path in index 0, 
	 * name in index 1 and possibly empty extension in index 2.
	 * 
	 * @return array with 3 possibly empty strings 
	 */
	public static String[] parseFilename(File file) {
		String path = file.getParent() != null ? file.getParent() : "";
		String wholeName = file.getName();
		String name, extension;
		if (wholeName.contains(".")) {
			name = wholeName.substring(0, wholeName.indexOf('.'));
			extension = wholeName.substring(wholeName.indexOf('.') + 1);
		} else {
			name = wholeName;
			extension = "";
		}
		
		return new String[] {path, name, extension};
	}
}
