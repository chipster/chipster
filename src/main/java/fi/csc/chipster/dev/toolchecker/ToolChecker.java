package fi.csc.chipster.dev.toolchecker;

import fi.csc.chipster.toolbox.SADLTool;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLParser.ParseException;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.io.File;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.nodes.Element;
import org.jsoup.select.Elements;


public class ToolChecker {
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
	 * @return Parameters found in the manual as a Set of Strings. Returns null if the file was not found.
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
			System.err.println("File not found: " + documentationPath);
			e.printStackTrace();
		}
		return null;
	}
	/**
	 * Updates manuals parameters according to new Parameters.
	 * @param newParameters, the string that is displayed on the manual.
	 * @param documentationPath, the path of the documentation
	 */
	private void writeParameters(Set<String> newParameters, String documentationPath) {
		File documentation = new File(documentationPath);
		Document doc;
		try {
			doc = Jsoup.parse(documentation, "UTF-8");
			Elements parameters = doc.select("h3:contains(Parameters) + ul");
			parameters.html("");
			// Modify parameters
			for (String p: newParameters) {
				parameters.prepend("<li>" + p + "</li>" );
			}
			// Write file
			FileWriter newFile = new FileWriter(documentation, false);
			newFile.write(doc.outerHtml());
			newFile.close();
			System.out.println(doc.outerHtml());
						
		} catch (IOException e) {
			e.printStackTrace();
		}
				
	}
		
	/**
	 * Returns null if there is a manual that has the same name as rScriptPath and manual has all the parameters, else returns an error message.
	 * If there is no manual, it prints the Script's name. Calls parametersError
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
	 * Assumes that the script and manual has the same name, except the file extension.
	 * Returns null if the parameters are valid. Calls parseR and parseHTML to obtain parameters.
	 * If parameters do not match returns an error string, white space errors are ignored
	 * @param filePath
	 * @param manualDir
	 * @return
	 */
	private String parametersErrors(Path filePath, String manualDir) {
		String fileName = filePath.getFileName().toString();
		fileName = fileName.substring(0,fileName.lastIndexOf('.'));		
		Set<String> manualParameters = parseHTML(manualDir + '/' + fileName + ".html");		
		Set<String> rParameters = parseR(filePath.toString());
		Set<String> filteredRParameters = new HashSet<String>();
		for (String p : rParameters){
			// Filter white space
			filteredRParameters.add(p.replaceAll("\\s+", " ").trim().replaceAll("\\.", ""));
		}
		
		if (manualParameters.isEmpty()) {
			if (rParameters.isEmpty()) {
				return null;
			} else {
				//writeParameters(filteredRParameters, manualDir + '/'+ fileName + ".html");
				return ("No parameters in manual" + "\n");
			}
		}
		else if (manualParameters.containsAll(filteredRParameters)) {
			return null;
		} else {
			//writeParameters(filteredRParameters, manualDir + '/'+ fileName + ".html");
			Set<String> rCopy = new HashSet<String>(filteredRParameters);
			System.err.println(filteredRParameters.toString());
			System.err.println(manualParameters.toString());
			filteredRParameters.removeAll(manualParameters);
			manualParameters.removeAll(rCopy);
			return ("Script's parameters mismatches: " + filteredRParameters.toString() + "\n" 
					+ "Manual's parameters mismatches: " + manualParameters.toString() + "\n");
		}
	}

	/**
	 * Calls manualError for each R script in the directory.
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
					//TODO replace hard coded
					errorMessage = manualError(child.toPath(), "/home/aoikari/git/chipster-tools/manual");
					if (errorMessage == null) {
		        		System.out.println("Manual is OK: " + child.toString());
		        		++ok;
					} else {
		        		System.out.println("MANUAL IS NOT OK: " + child.toString() + "\n" + errorMessage );
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
		// TODO: remove hard coded paths
		new ToolChecker().scanDir("/home/aoikari/git/chipster-tools/tools/ngs/R");		
	}
}