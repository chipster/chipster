package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.HashMap;

public class GeneSet extends HashMap<String, Gene>{
	
	public void addExon(Exon exon, String geneId, String transcId,
			String geneName, String transcName, String biotype) {
		
		Gene gene;
		
		if ((gene = this.get(geneId)) == null) {
			gene = new Gene(geneName, biotype, geneId);
			this.put(geneId, gene);
		}
		
		gene.addExon(exon, geneId, transcId, transcName);
	}
	
	public Gene getGene(String name) {
		
		for (Gene gene : this.values()) {
			if (name.equalsIgnoreCase(gene.getName())) {
				return gene;
			}
		}
		return null;
	}
}
