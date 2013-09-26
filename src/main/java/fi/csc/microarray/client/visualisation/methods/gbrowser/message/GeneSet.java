package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.HashMap;
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
	 * Iterate over exons and add them to this GeneSet. Iterator.remove() is called
	 * for all exons that don't intercept with the filterRegion.
	 * 
	 * @param iterator
	 * @param filterRegion
	 */
	public void add(Iterator<Exon> iterator, Region filterRegion) {	
		
		while (iterator.hasNext()) {
			
			Exon exon = iterator.next();
			
			if (!filterRegion.intersects(exon.getRegion())) {
				iterator.remove();
			}
			
			addExon(exon, exon.getGeneId(), exon.getTranscriptId(), exon.getGeneName(), exon.getTranscName(), exon.getBiotype());			
		}
	}
}
