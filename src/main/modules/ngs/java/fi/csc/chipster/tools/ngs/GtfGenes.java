package fi.csc.chipster.tools.ngs;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.GtfToFeatureConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneSet;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.GtfLineParser;

public class GtfGenes {
	public static void main(String[] args) throws Exception {

		try {
		
		File in = new File(args[0]);
		File ids = new File(args[1]);
		File names = new File(args[2]);
		
		System.out.println("Parsing genes from GTF file " + in + "...");
		
		ids.delete();
		names.delete();
				
		getGenes(in, ids, names);
		
		System.out.println("DONE");
		
		} catch (Exception e) {
			e.printStackTrace();
						
			System.out.println(
					"usage: \n" +
					"  GtfGenes input-gtf gene-id-output gene-name-output\n" +
					"example:\n " +
					"  java -cp chipster-3.0.0.jar fi.csc.chipster.tools.ngs.GtfGenes Homo_sapiens.GRCh37.70.gtf Homo_sapiens.GRCh37.70.gene-id.search.tsv Homo_sapiens.GRCh37.70.gene-name.search.tsv");
		}				
	}

	private static void getGenes(File in, File ids, File names) {
		GtfLineParser parser = new GtfLineParser();
		LinkedList<Exon> exons = new LinkedList<Exon>();
		
		String lastChr = null;		
		
		/* Reading whole file at once consumes too much memory, so write out results
		 * always on chromosome change. Cutting anywhere else would be more difficult, 
		 * because gene must not overlap the cut position. 
		 */
		
		try (BufferedReader br = new BufferedReader(new FileReader(in))) {
			String line;
			while ((line = br.readLine()) != null) {
				Exon exon = GtfToFeatureConversion.parseLine(parser, line);
				
				if (exon == null) {					
					// all unrecognized feature types are skipped (e.g. start_codon and stop_codon)
					continue;
				}
				
				String currentChr = exon.getRegion().start.chr.toString();
				if (lastChr != null && !currentChr.equals(lastChr)) {
					write(ids, names, exons);
					exons.clear();					
				}
				lastChr = currentChr;
				
				exons.add(exon);			
			}
			// write out last chromosome
			write(ids, names, exons);
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private static void write(File ids, File names, LinkedList<Exon> exons) {		

		GeneSet geneSet = new GeneSet();				
		geneSet.add(exons.iterator(), null);
		
		String chr = null;
						
		try ( BufferedWriter idWriter = new BufferedWriter(new FileWriter(ids, true));
				BufferedWriter nameWriter = new BufferedWriter(new FileWriter(names, true))) {

			for (Gene gene : geneSet.values()) {
			
				String name = gene.getName();
				String id = gene.getId();
				
				chr = gene.getRegion().start.chr.toString();
				
				if (id != null) {
					idWriter.write(
							chr + "\t" + 
									gene.getRegion().start.bp + "\t" +  
									gene.getRegion().end.bp + "\t" +  
									id + "\n");
				}
				
				if (name != null) {
					nameWriter.write(
							chr + "\t" + 
									gene.getRegion().start.bp + "\t" +  
									gene.getRegion().end.bp + "\t" + 
									name + "\n");
				}
			}
			
			// does try-with-resources do this already?
			idWriter.flush();
			nameWriter.flush();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}		
		System.out.println("chromosome " + chr);
	}
}
