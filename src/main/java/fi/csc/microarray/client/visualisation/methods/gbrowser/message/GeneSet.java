package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

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
	
	/**
	 * Add all exons in the given set and remove exons that don't intercept with the fitlerRegion.
	 * 
	 * @param exons
	 * @param filterRegion
	 */
	public void add(HashSet<Exon> exons, Region filterRegion) {
		
		Iterator<Exon> exonIter = exons.iterator();
		
		while (exonIter.hasNext()) {
			
			Exon exon = exonIter.next();
			
			if (!filterRegion.intersects(exon.getRegion())) {
				exonIter.remove();
			}
			
			addExon(exon, exon.getGeneId(), exon.getTranscriptId(), exon.getGeneName(), exon.getTranscName(), exon.getBiotype());
		}
	}
}
