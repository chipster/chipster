package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Collection;
import java.util.HashMap;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class GeneSet extends TreeMap<Region, Gene>{
	
	//Genes mapped according to their id's to find them easily when the data is constructed
	HashMap<String, Gene> geneIdMap = new HashMap<String, Gene>();
	
	HashMap<String, Transcript> transcriptIdMap = new HashMap<String, Transcript>();
	
	private Gene gene;
	private Transcript transcript;
	
	public void addExon(Exon exon, String geneId, String transcId,
			String geneName, String transcName, String biotype) {
		
		if ((gene = geneIdMap.get(geneId)) == null) {
			gene = new Gene(geneName, biotype);
			geneIdMap.put(geneId, gene);;
		}
		
		if((transcript = transcriptIdMap.get(transcId)) == null) {
			transcript = new Transcript(transcName);
			transcriptIdMap.put(transcId, transcript);
		}
		
		transcript.addExon(exon);
		
		gene.addTranscript(transcript);
	}
	
	public void prepareForReading() {
		
		for (Gene gene : geneIdMap.values()) {
			this.put(gene.getRegion(), gene);
			
			gene.prepareForReading();
		}
						
		geneIdMap = null;
		transcriptIdMap = null;
	}

	public Collection<Gene> getGenes(Region region) {
		
		return this.subMap(new Region(region.start, region.start), new Region(region.end, region.end)).values();
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
