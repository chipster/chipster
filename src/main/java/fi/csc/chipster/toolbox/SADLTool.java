package fi.csc.chipster.toolbox;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintWriter;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;
import fi.csc.microarray.util.IOUtils;

/**
 * Tool for SADL related tasks.
 * 
 * @author Aleksi Kallio
 */
public class SADLTool {

	public static class ParsedScript {
		public String SADL = "";
		public String source = "";
		public String code = "";
	}
	
	private String comment;
	
	public SADLTool(String comment) {
		this.comment = comment;
	}


	public ParsedScript parseScript(InputStream rScriptSource) throws IOException {

		ParsedScript parsedScript = new ParsedScript();		
		BufferedReader in = null;
		try {
			in = new BufferedReader(new InputStreamReader(rScriptSource));
			boolean inHeaderCommentBlock = true;
			for (String line = in.readLine(); line != null; line = in.readLine()) {
				if (!line.startsWith(comment)) {
					inHeaderCommentBlock = false;
				}
				if (inHeaderCommentBlock) {
					String strippedLine = line.substring(comment.length());
					if (strippedLine.startsWith(" ")) {
						strippedLine = strippedLine.substring(1); // remove space after comment symbol
					}
					parsedScript.SADL += strippedLine + "\n";
					
				} else {
					parsedScript.code += line + "\n";
				}
				
				parsedScript.source += line + "\n";
			}
		} finally {
			IOUtils.closeIfPossible(in);
		}
		
		return parsedScript;
	}
	
	public String toScriptString(ParsedScript parsedScript) {
		String string = "";
		
		// Do SADL part
		for (String sadlRow : parsedScript.SADL.split("\n")) {
			string += comment + " " + sadlRow + "\n";
		}
		
		// Do rest
		string += parsedScript.code;
		
		return string;
	}
	
	public void convertToSADL(File scriptFile) throws IOException, ParseException {
		FileInputStream in = null;
		ParsedScript parsedScript;
		try {
			// Load the script
			in = new FileInputStream(scriptFile);
			parsedScript = parseScript(in);
			
			// Parse it
			SADLDescription sadl = new ChipsterSADLParser().parse(parsedScript.SADL);
			
			// Do necessary tweaks to meta inputs/outputs
			for (Input input : sadl.getInputs()) {
				if ("phenodata.tsv".equals(input.getName().getID())) {
					input.setMeta(true);
				}
			}
			for (Output output : sadl.getOutputs()) {
				if ("phenodata.tsv".equals(output.getName().getID())) {
					output.setMeta(true);
				}
			}

			// Generate SADL from parsed content
			parsedScript.SADL = sadl.toString();
			
		} finally {
			IOUtils.closeIfPossible(in);
		}
		
		// Rewrite the file
		PrintWriter out = null;
		try {
			out = new PrintWriter(new FileOutputStream(scriptFile));
			out.print(toScriptString(parsedScript));
			
		} finally {
			IOUtils.closeIfPossible(out);
		}
	}
	
	public static void main(String[] args) throws Exception {
		DirectoryLayout.initialiseUnitTestLayout();
		convertBeanshell();
	}


	public static void convertACGH() throws Exception {
		for (String moduleDir : new String[] { "src/main/modules/microarray/R-2.12", "src/main/modules/ngs/R-2.12" }) {
			for (File file : new File(moduleDir).listFiles()) {
				if (file.getName().endsWith(".R")) {
					System.out.println(file.getCanonicalPath());
					new SADLTool("#").convertToSADL(file);
				}
			}
		}
	}
	
	public static void convertSomeScripts() throws Exception {
		File files[] = new File[] {
			new File("src/main/modules/microarray/R-2.12/annotate-agilent-miRNA.R"),
			new File("src/main/modules/ngs/R-2.12/ngs-annotate-miRNA-targets.R"),
			new File("src/main/modules/ngs/R-2.12/ngs-correlation-analysis-mirna.R"),
			new File("src/main/modules/ngs/R-2.12/ngs-extract-target-gene-expression.R"),
			new File("src/main/modules/ngs/R-2.12/ngs-filter-annotations.R"),
			new File("src/main/modules/ngs/R-2.12/ngs-filter-results-column.R"),
			new File("src/main/modules/ngs/R-2.12/ngs-pathways-mirna-hyperg-go.R"),
			new File("src/main/modules/ngs/R-2.12/ngs-pathways-mirna-hyperg-kegg.R"),
			new File("src/main/modules/ngs/R-2.12/ngs-up-down-analysis-mirna.R")
		};
		
		for (File file : files) {
			if (file.getName().endsWith(".R")) {
				System.out.println(file.getCanonicalPath());
				new SADLTool("#").convertToSADL(file);
			}
		}
	}

	public static void convertBeanshell() throws Exception {
		for (File file : new File("src/main/modules/microarray/bsh").listFiles()) {
			if (file.getName().endsWith(".bsh")) {
				System.out.println(file.getCanonicalPath());
				new SADLTool("//").convertToSADL(file);
			}
		}
	}


}
