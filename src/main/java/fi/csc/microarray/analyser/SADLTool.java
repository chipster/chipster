package fi.csc.microarray.analyser;

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
			SADLDescription sadl = new ChipsterSADLParser().parse(parsedScript.SADL, scriptFile.getName());
			
			// Do necessary tweaks to meta inputs/outputs
			for (Input input : sadl.inputs()) {
				if ("phenodata.tsv".equals(input.getName().getID())) {
					input.setMeta(true);
				}
			}
			for (Output output : sadl.outputs()) {
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
		convertACGH();
	}


	private static void convertACGH() throws Exception {
//		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/add-cytobands.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/calculate-aberration-frequencies.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/cluster-acgh.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/convert-cn-probes-to-genes.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/count-overlapping-cnvs.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/detect-common-copy-number-aberration-regions.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/detect-copy-number-aberrations.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/fetch-probe-positions-from-cangem.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/import-from-cangem.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/match-cn-and-expression-probes.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/merge-datasets.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/pathways-acgh-hyperg-go.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/plot-cgh-profile.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/plot-cn-induced-expression-profile.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/plot-cn-induced-gene-expression.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/sample-size-with-bh.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/search-coordinates.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/smooth-acgh.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/stat-acgh.R"));
		new SADLTool("#").convertToSADL(new File("src/main/tools/R-2.12/test-for-cn-induced-differential-expression.R"));
	}

}
