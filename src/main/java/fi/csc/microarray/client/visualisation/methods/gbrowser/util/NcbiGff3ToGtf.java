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
 * This implementation is not tested, because a prerelease of the annotation was found from Ensembl.
 * 
 * @author klemela
 */
public class NcbiGff3ToGtf {

	public static void main(String[] args) throws IOException {

		File in = new File("/home/klemela/chipster/gb-species/chicken/ref_Gallus_gallus-4.0_top_level.gff3");
		File out = new File("/home/klemela/chipster/gb-species/chicken/ref_Gallus_gallus-4.0_top_level.gtf");

		FileInputStream inStream = new FileInputStream(in);
		BufferedReader reader = new BufferedReader(new InputStreamReader(inStream));
		String line = null;

		BufferedWriter writer = new BufferedWriter(new FileWriter(out));		

		String humanReadableChr = null;
		String idOfHumanReadableChr = null;

		int exonCounter = -1;
		String transcriptIdOfExonCounter = null;

		int cdsCounter = -1;
		String transcriptIdOfCdsCounter = null;

		try {

			while ((line = reader.readLine()) != null)   {

				if (line.startsWith("#")) {
					continue;
				}

				String[] cols;

				Map<String, String> ids;

				//Read input
				cols = line.split("\t");

				String chr = cols[0];
				String biotype = cols[1];
				String feature = cols[2];
				String exonStart = cols[3];
				String exonEnd = cols[4];
				String frame = cols[5];
				String strand = cols[6];
				String score = cols[7];


				ids = parseIds(cols[8]);

				//Chromosome 
				if ("region".equals(feature)) {
					idOfHumanReadableChr = chr;
					humanReadableChr = ids.get("chromosome");
					continue;
				}
				
				//What is this?
				if ("match".equals(feature)) {					
					continue;
				}

				if (!chr.equals(idOfHumanReadableChr)) {
					System.err.println("Unexpected chromosome identifier on line: " + line);
					break;
				}
				chr = humanReadableChr;

				String geneId = ids.get("Dbxref");
				ids.remove("Dbxref");

				String outLine = null;

				if ("exon".equals(feature)) {

					String geneName = ids.get("gene");
					ids.remove("gene");			

					String exonId = ids.get("ID");
					ids.remove("ID");			
					String transcId = ids.get("transcript_id");
					ids.remove("transcript_id");
					
					if (transcId == null) {
						transcId = exonId;
					}

					//Convert some fields
					biotype = "biotype";

					if (transcId.equals(transcriptIdOfExonCounter)) {
						exonCounter++;
					} else {
						exonCounter = 0;
						transcriptIdOfExonCounter = transcId;
					}

					String exonIndex = "" + exonCounter;					


					outLine = 
							chr + "\t" + 
									biotype + "\t" + 
									feature + "\t" + 
									exonStart + "\t" + 
									exonEnd + "\t" + 
									score + "\t" + 
									strand + "\t" +
									frame + "\t" +
									" gene_id \"" + geneId + "\";" + 
									" transcript_id \"" + transcId + "\";" +
									" exon_number \"" + exonIndex + "\";" +
									" gene_name \"" + geneName + "\";" + 
									" exonId \"" + exonId + "\";";


				} else if ("CDS".equals(feature)) {
			
					String transcId = ids.get("protein_id");
					ids.remove("protein_id");
					
					if (transcId == null) {
						transcId = ids.get("ID");
					}


					//Convert some fields
					biotype = "biotype";

					if (transcId.equals(transcriptIdOfCdsCounter)) {
						cdsCounter++;
					} else {
						cdsCounter = 0;
						transcriptIdOfCdsCounter = transcId;
					}

					String exonIndex = "" + cdsCounter;					

					outLine = 
							chr + "\t" + 
									biotype + "\t" + 
									feature + "\t" + 
									exonStart + "\t" + 
									exonEnd + "\t" + 
									score + "\t" + 
									strand + "\t" +
									frame + "\t" +
									" gene_id \"" + geneId + "\";" + 
									" transcript_id \"" + transcId + "\";" +
									" exon_number \"" + exonIndex + "\";";

				}


				for (Map.Entry<String, String> entry : ids.entrySet()) {
					outLine += " " + entry.getKey() + " \"" + entry.getValue() + "\";";
				}

				writer.write(outLine + "\n");
			}

		} catch (Exception e) {
			e.printStackTrace();
			System.out.println(line);
		}

		inStream.close();
		writer.flush();
		writer.close();
	}

	public static Map<String, String> parseIds(String ids) {		

		String[] split = ids.split(";");
		Map<String, String> result = new HashMap<String, String>();

		String key = null;
		String value = null;
		int indexOfQuotationMark = 0;

		for (int i = 0; i < split.length; i++) {

			indexOfQuotationMark = split[i].indexOf("=");

			key = split[i].substring(0, indexOfQuotationMark);
			value = split[i].substring(indexOfQuotationMark + 1, split[i]
					.length());

			result.put(key, value);
		}

		return result;
	}
}
