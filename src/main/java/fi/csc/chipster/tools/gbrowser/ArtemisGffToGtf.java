package fi.csc.chipster.tools.gbrowser;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

public class ArtemisGffToGtf {

	public static void main(String[] args) throws IOException {
		
//		File in = new File("/home/klemela/.chipster/yersinia/N916Ysi.gff");
//		File out = new File("/home/klemela/.chipster/yersinia/N916Ysi.gtf");
		
		File in = new File("/home/klemela/.chipster/yersinia/R1-RT.gff");
		File out = new File("/home/klemela/.chipster/yersinia/R1-RT.gtf");

		FileInputStream inStream = new FileInputStream(in);
		BufferedReader reader = new BufferedReader(new InputStreamReader(inStream));
		String line;

		BufferedWriter writer = new BufferedWriter(new FileWriter(out));		

		while ((line = reader.readLine()) != null)   {

			if ("##FASTA".equals(line)) {
				break;
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
			
			//R1-RT file has a duplicate gene row for every exon
			if ("gene".equals(feature)) {
				continue;
			}			

			ids = parseIds(cols[8]);

			String exonId = ids.get("ID");
			ids.remove("ID");			
			String geneName = ids.get("gene");
			ids.remove("gene");


			//Convert some fields
			chr = "Chromosome";
			biotype = "biotype";
			String geneId = exonId;
			String transcId = exonId;
			String exonIndex = "1";

			String outLine = 
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

			if (geneName != null) {
				outLine += "gene_name \"" + geneName + "\"; ";
			}


			for (Map.Entry<String, String> entry : ids.entrySet()) {
				outLine += " " + entry.getKey() + " \"" + entry.getValue() + "\";";
			}

			writer.write(outLine + "\n");
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
