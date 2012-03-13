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
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

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
	 * Lists all files (for which isDirectory() returns false) 
	 * under this file or directory and its subdirectories.
	 * If called with a non-directory input, the input itself
	 * is returned.
	 * 
	 * @param file directory (or file) to be recursed
	 * 
	 * @return list of Files
	 */
	public static List<File> listFilesRecursively(File file) {
		LinkedList<File> files = new LinkedList<File>();
		
		if (file.isDirectory()) {
			// dir, recurse into it and combine result lists
			for (File subFile : file.listFiles()) {
				files.addAll(listFilesRecursively(subFile));
			}
			
		} else {
			// file, add and return this
			files.add(file);
		}
		
		return files;		
	}
	
	
	/**
	 * Deletes a file or a directory recursively. Deletes directory links, does not go 
	 * into them recursively.
	 * 
	 * @param dir directory or file to be deleted
	 * @return true if deleting was successful, false file does not exist or deleting it failed
	 * @throws IOException 
	 */
	public static boolean delTree(File dir) throws IOException {

		// Just try to delete the file first
		// Will work for normal files, empty dirs and links (dir or file)
		// Avoids need for dealing with links later on
		if (dir.delete()) {
			return true;
		} 
		
		// Directory
		else if (dir.isDirectory()) {
			for (File file : dir.listFiles()) {
				delTree(file);
			}

			// Dir should be empty now
			return dir.delete();
		} 
		
		// Could not delete, not a directory, no can do
		else {
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

	public static boolean sweepAndCreateDirectory(File dir) throws IOException {
		
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
	 * Walks the baseDir recursively and deletes files and directories older than cutoff. 
	 * When traversing directories, does not follow symbolic links.
	 * 
	 * If a directory is old but contains files (which are not too old), it is not deleted.
	 * 
	 * TODO better problem handling?
	 * 
	 * @param baseDir
	 * @param cutoff milliseconds 
	 * @throws IOException 
	 */
	public static void cleanOldFiles(File baseDir, long cutoff ) throws IOException {
		walkAndDelete(baseDir, new AgeFileFilter(System.currentTimeMillis() - cutoff));
	}

	/**
	 * Walks the baseDir recursively and deletes files that match filter. When traversing directories, does not follow symbolic links.
	 * 
	 * @param baseDir
	 * @param filter
	 * @throws IOException 
	 */
	public static void walkAndDelete(File baseDir, IOFileFilter filter) throws IOException {
		File[] files = baseDir.listFiles((FileFilter)new OrFileFilter(DirectoryFileFilter.INSTANCE, filter));
		
		if (files == null) {
			return;
		}
		
		for (File f : files) {
			
			// Just try to delete it first
			// Will work for normal files, empty dirs and links (dir or file)
			// Avoid need for dealing with links later on
			if (f.delete()) {
				continue;
			} else if (f.isDirectory()) {
				
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

	/**
	 * Converts File to URL, but does not throw checked MalformedURLException.
	 * 
	 * @return URL or if it is malformed, then null
	 */
	public static URL toUrl(File file) {
		try {
			return file.toURI().toURL();
		} catch (MalformedURLException e) {
			return null;
		}
	}
	
	/**
	 * Find files in a given directory whose filenames match given regex.
	 */
	public static File[] findFiles(File dir, String regex) {
	    
	    class RegexFileFilter implements FilenameFilter {
	        private String regex;
	        
	        public RegexFileFilter(String regex) {
                this.regex = regex;
            }

            @Override
            public boolean accept(File dir, String name) {
                return name.matches(regex);
            }
	        
	    }
	    
	    return dir.listFiles(new RegexFileFilter(regex));
	}
	
	
	/**
	 * Tries to create a symbolic link. Returns true is link was created and false otherwise (links not supported 
	 * on the platform, IO error, ...). Canonical paths of from and to Files are used when creating the link.
	 *  
	 * @param from
	 * @param to
	 * @return true iff link created successfully
	 */
	public static boolean createSymbolicLink(File from, File to) {
		Process process = null;
		try {
			process = Runtime.getRuntime().exec( new String[] { "ln", "-s", from.getCanonicalPath(), to.getCanonicalPath() } );
			process.waitFor();
			return true;
			
		} catch (Exception e) {
			return false;
			
		} finally {
			if (process != null) {
				process.destroy();
			}
		}
	}
	
	public static void main(String[] args) throws IOException {
		new File("/home/akallio/link_session").delete();
		new File("/home/akallio/link_test").delete();
	}
}
