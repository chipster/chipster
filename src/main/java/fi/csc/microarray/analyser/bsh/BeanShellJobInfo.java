package fi.csc.microarray.analyser.bsh;

import java.io.File;
import java.util.HashMap;

/**
 * Wrapper for the job information passed to the BeanShell scripts.
 *
 * 
 * @author hupponen
 *
 */
public class BeanShellJobInfo {

	public File workDir;
	
	public HashMap<String, String> parameters = new HashMap<String,String>();
	
}
