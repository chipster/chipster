package fi.csc.microarray.util;

import java.io.File;

import javax.swing.filechooser.FileFilter;

/**
 * A file filter used by different file choosers of this application.
 * Can be instantiated with any list of allowed file extensions (given
 * as a String array) to restrict what kind of files are shown, for
 * example, when selecting a file for importing microarray data.
 * 
 * @author Janne KÃ¤ki
 *
 */
public class GeneralFileFilter extends FileFilter {

	private String description;
	private String[] extensions;
	
	/**
	 * Creates a new file filter.
	 * 
	 * @param description A short description telling what kind of files
	 * 					  this filter is intended to accept. Is shown in
	 * 					  the filter combobox of the file chooser.
	 * @param extensions An array of all file extensions (without the
	 *                   separator period) that this filter should accept.
	 */
	public GeneralFileFilter(String description, String[] extensions) {
		this.description = description;
		this.extensions = extensions;
	}
	
	/**
	 * A method for determining whether this filter should accept the
	 * given file or not (accepting meaning that it will be displayed in
	 * the file chooser that has this filter assigned). A filter will
	 * accept all directories plus all files with a "known" extension
	 * (one given to its constructor).
	 * 
	 * @param file The file to be evaluated.
	 * @return True if the file is accepted, false if not.
	 */
    public boolean accept(File file) {
		if (file.isDirectory()) {
		    return true;
		}
		String extension = getExtension(file);
		if (extension != null) {
			for (String s : extensions) {
			    if (s.equals(extension)) {
			    	return true;
			    }
			}
		}
		return false;
    }

    /**
     * A helper method from extracting the extension (the substring following
     * the period, '.', in the file name) from a file.
     * 
     * @param file The file whose extension might be interesting.
     * @return The extension as a String (if any can be found).
     */
    private String getExtension(File file) {
		String extension = null;
		String filename = file.getName();
		int i = filename.lastIndexOf('.');
		if ( i > 0 && i < filename.length()-1 ) {
		    extension = filename.substring(i+1).toLowerCase();
		}
		return extension;
    }
	    
    /**
     * @return The description of this filter (one given as a parameter to
     * 		   the constructor).
     */
    public String getDescription() {
    	return description;
    }
}
