package fi.csc.chipster.dev.toolchecker;

import fi.csc.chipster.toolbox.SADLTool;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLParser.ParseException;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FilenameFilter;
import java.io.IOException;
import java.nio.file.FileSystems;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.PathMatcher;
import java.nio.file.Paths;
import java.io.File;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;

import com.fasterxml.jackson.databind.ObjectMapper;

public class ToolChecker {
	private static String errorMessage;


	/**
	 * Parses the arguments' names of a R script.
	 * @param rPath, path to R script.
	 * @return The parameters of the script in a String Set.
	 */
	private Set<String> parseR(String rPath) {		
			File file = new File(rPath);
			FileInputStream fileStream;
			Set<String> parametersList = new HashSet<String>();
			try {
				fileStream = new FileInputStream(file);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
				return parametersList;
			}
			try {
				SADLTool.ParsedScript script = new SADLTool("#").parseScript(fileStream);
				SADLParser parser = new SADLParser();
				SADLDescription description = parser.parse(script.SADL);
				List<Parameter> parameters = description.getParameters();
				for (Parameter p : parameters) {
					parametersList.add(p.getName().getDisplayName());
				}
			} catch (ParseException | IOException e) {
				System.err.println(e);
			}
			return parametersList; 		
	}
	
	/**
	 * Parses the arguments's names from a given manual page.
	 * @param documenationPath, path of the manual.
	 * @return Parameters found in the manual in a String Set.
	 */
	private Set<String> parseHTML(String documentationPath) {
		File documentation = new File(documentationPath);
		Document doc;
		try {
			doc = Jsoup.parse(documentation, "UTF-8");
			// Finds list items from a list that is right after a h3 header that contains "Parameters"
			Elements parameters = doc.select("h3:contains(Parameters) + ul > li");
			Set<String> parametersList = new HashSet<String>();
			for ( Element p :parameters) {
				// Split the text from (, [, . ,and : before those chars is the name that we want
				String[] tokens = p.text().split("\\(|\\[|\\.|\\:");
				parametersList.add(tokens[0].trim());
			}
			return parametersList;			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.err.println("File not found: " + documentationPath);
			e.printStackTrace();
		}
		return null;
	}
		
	/**
	 * Returns null if there is a manual that has the same name as rScriptPath and manual has all the parameters.
	 * If there is no manual, it prints the Script's name.
	 * @param rScriptPath, path of R script.
	 * @param manualDir, manual directory.
	 * @return
	 */
	private String manualError(Path rScriptPath, String manualDir) {		
		String fileName = rScriptPath.getFileName().toString();
		fileName = fileName.substring(0,fileName.lastIndexOf('.'));
		String filePath = manualDir + '/' + fileName + ".html";
		if (new File(filePath).isFile()) {
			return parametersErrors(rScriptPath, manualDir);
		} else {
			return ("MANUAL NOT FOUND: " + fileName);
		}
	}
	 
	/**
	 * Assumes that the script and manual has the same name, except file extension.
	 * Returns null if the parameters are valid.
	 * If parameters do not match returns a error string
	 * @param filePath
	 * @param manualDir
	 * @return
	 */
	private String parametersErrors(Path filePath, String manualDir) {
		String fileName = filePath.getFileName().toString();
		fileName = fileName.substring(0,fileName.lastIndexOf('.'));		
		Set<String> manualParameters = parseHTML(manualDir + '/' + fileName + ".html");		
		Set<String> rParameters = parseR(filePath.toString());
		if (manualParameters.isEmpty()) {
			if (rParameters.isEmpty()) {
				return null;
			} else {
				return ("No parameters in manual" + "\n");
			}
		}
		else if (manualParameters.containsAll(rParameters)) {
			return null;
		} else {
			Set<String> rCopy = new HashSet<String>(rParameters);
			rParameters.removeAll(manualParameters);
			manualParameters.removeAll(rCopy);
			return ("Script's parameters mismatches: " + rParameters.toString() + "\n" 
					+ "Manual's parameters mismatches: " + manualParameters.toString() + "\n");
		}
	}

	/**
	 * Calls search manual and verifyParameters for each R script in the directory.
	 * Prints out is the manual OK or not.
	 * @param rDirPath, path of the R script directory
	 */
	private void scanDir(String rDirPath) {		
		int ok = 0;
		int error = 0;
		File dir = new File(rDirPath);
		File[] directoryListing = dir.listFiles();
		String errorMessage = null;
		if(directoryListing != null) {
			for (File child : directoryListing) {
				if (child.isFile()) {
					errorMessage = manualError(child.toPath(), "/home/aoikari/git/chipster-tools/manual");
					if (errorMessage == null) {
		        		System.out.println("Manual is OK: " + child.toString());
		        		++ok;
					} else {
		        		System.err.println("MANUAL IS NOT OK: " + child.toString() + "\n" + errorMessage );
		        		++error;
		        	}					
				}				
			}
		} else {
			System.err.println("No directory: " + rDirPath);
		}		
		System.out.println("Manuals OK: " + ok + " Manuals not OK: " + error);		
	}
		
	public static void main(String[] args) {		
		new ToolChecker().scanDir("/home/aoikari/git/chipster-tools/tools/ngs/R");		
	}
}