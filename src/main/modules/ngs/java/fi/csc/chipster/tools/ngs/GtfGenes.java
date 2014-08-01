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
		File out = new File(args[1]);
		
		System.out.println("Parsing genes from GTF file " + in + " to " + out + "...");
		
		out.delete();
				
		getGenes(in, out);
		
		System.out.println("DONE");
		
		} catch (Exception e) {
			e.printStackTrace();
						
			System.out.println(
					"usage: \n" +
					"  GtfGenes <file-in> <file-out>\n" +
					"example:\n " +
					"  java -cp chipster-3.0.0.jar fi.csc.chipster.tools.ngs.GtfGenes Homo_sapiens.GRCh37.70.gtf Homo_sapiens.GRCh37.70.gene.tsv");
		}				
	}

	private static void getGenes(File in, File out) {
		GtfLineParser parser = new GtfLineParser();
		LinkedList<Exon> exons = new LinkedList<Exon>();
		
		String lastChr = null;			
		
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
					write(out, exons);
					exons.clear();					
					System.out.println("chromosome " + lastChr);
				}
				lastChr = currentChr;
				
				exons.add(exon);			
			}
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}
	}

	private static void write(File out, LinkedList<Exon> exons) {		

		GeneSet geneSet = new GeneSet();				
		geneSet.add(exons.iterator(), null);		
						
		try ( BufferedWriter writer = new BufferedWriter(new FileWriter(out, true))) {

			for (Gene gene : geneSet.values()) {
			
				String name = gene.getName();
				String id = gene.getId();
				
				name = name == null ? "" : name;
				id = id == null ? "" : id;
				
				writer.write(
						gene.getRegion().start.chr + "\t" + 
						gene.getRegion().start.bp + "\t" +  
						gene.getRegion().end.bp + "\t" + 
						name + "\t" + 
						id + "\n");
			}
			writer.flush();
			
		} catch (IOException e) {
			e.printStackTrace();
			System.exit(1);
		}		
	}
}
