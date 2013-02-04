package fi.csc.microarray.client.visualisation.methods.gbrowser.util;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

/**
 * @author klemela
 */
public class GenbankToGtf {

	private static final String START = "start";
	private static final String END = "end";
	private static final String STRAND = "strand";

	public static void main(String[] args) throws IOException {

		File in = new File("/home/klemela/chipster/gb-species/human-mitochondrion/sequence.gb");
		File out = new File("/home/klemela/chipster/gb-species/human-mitochondrion/sequence.gtf");

		String chr = "MT";

		FileInputStream inStream = new FileInputStream(in);
		BufferedReader reader = new BufferedReader(new InputStreamReader(inStream));
		String line = null;

		BufferedWriter writer = new BufferedWriter(new FileWriter(out));		

		try {

			while ((line = reader.readLine()) != null)   {
				if (line.startsWith("FEATURES")) {
					break;
				}
			}

			String col1 = null;
			String col2 = null;
			String feature = "";
			Map<String, String> ids = new HashMap<String, String>();
			String featureQualifier = "";

			while ((line = reader.readLine()) != null)   {

				if (line.startsWith("ORIGIN")) {
					//Rest of the file is sequences, conversion is ready
					parseFeatureQualifier(featureQualifier, ids);
					writeLine(ids, feature, chr, writer);
					break;
				}

				if (line.length() <= 20) {
					System.err.println("Skipping too short line: \"" + line + "\"");
					continue;					
				}
				
				col1 = line.substring(0, 20).trim();
				col2 = line.substring(20, line.length()).trim();


				if (col1.length() > 0) {
					
					//The line has a new feature, save previous feature first (unless this is a first one)
					parseFeatureQualifier(featureQualifier, ids);
					writeLine(ids, feature, chr, writer);
					
					feature = parseFeatureHeader(col1, col2, ids);

					continue;
				} else {
					if (col2.startsWith("/")) {	
						
						parseFeatureQualifier(featureQualifier, ids);
						featureQualifier = col2;
					} else {
						featureQualifier += " " + col2;
					}
				}
			}
			//Write last line
			parseFeatureQualifier(featureQualifier, ids);
			writeLine(ids, feature, chr, writer);

		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(line);
		}

		inStream.close();
		writer.flush();
		writer.close();
	}

	private static void parseFeatureQualifier(String line, Map<String, String> ids) {
		if (line.length() > 0) {
			String[] keyAndValue = line.substring(1, line.length()).split("=");

			String key = keyAndValue[0];
			String value = keyAndValue[1];

			ids.put(key, value);
		}
	}

	private static String parseFeatureHeader(String column1, String column2, Map<String, String> ids) {

			String feature = column1;

			if (column2.contains("complement")) {
				ids.put(STRAND, "-");
				column2 = column1.replace("complement(", "").replace(")", "");
			} else {
				ids.put(STRAND, "+");
			}
			
			if (column2.contains("..")) {
				//Region
				String[] region = column2.split("\\.\\.");
				String start = region[0].replace("<", "");
				String end = region[1].replace(">", "");					
				ids.put(START, start);
				ids.put(END, end);
			} else {
				//Just single coordinate
				ids.put(START, column2);
				ids.put(END, column2);
			}
			
			return feature;
	}

	private static void writeLine(Map<String, String> ids, String feature, String chr, BufferedWriter writer) throws IOException {
		if (ids.size() > 0) {
			String biotype = "biotype";
			feature = "GenBank " + feature;
			String start = ids.get(START);
			ids.remove(START);
			String end = ids.get(END);
			ids.remove(END);
			String score = ".";
			String strand = ids.get(STRAND);
			ids.remove(STRAND);

			//How to calculate this?
			String frame = ".";

			String geneId = ids.get("gene");
			String transcId = geneId;

			String outLine = 
					chr + "\t" + 
							biotype + "\t" + 
							feature + "\t" + 
							start + "\t" + 
							end + "\t" + 
							score + "\t" + 
							strand + "\t" +
							frame + "\t";
			
			if (geneId != null) {
				outLine += " gene_id \"" + geneId + "\";";
			}
			
			if (transcId != null) {
				outLine += " transcript_id \"" + transcId + "\";";
			}

			outLine += mapToString(ids);

			System.out.println(outLine);
			
			writer.write(outLine + "\n");

			ids.clear();
		}
	}

	public static String mapToString(Map<String, String> map) {
		String string = "";
		for (Map.Entry<String, String> entry : map.entrySet()) {
			string += " " + entry.getKey() + " \"" + entry.getValue() + "\";";
		}	
		return string;
	}
}
