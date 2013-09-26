package fi.csc.microarray.client.visualisation.methods.gbrowser.message;


public class Exon implements Comparable<Exon> {
	
	private Region region;
	private Feature feature;
	private int exonNumber;
	private Transcript transcript;

	private String geneId;
	private String transcriptId;
	private String geneName;
	private String transcName;
	private String biotype;

	public static enum Feature {

		UNRECOGNIZED(null),
		TRANSCRIPT("transcript"), //Cufflinks gtf
		EXON("exon"),
		CDS("CDS"), 
		START_CODON("start_codon"),
		STOP_CODON("stop_codon");

		private String id;

		Feature(String id) {
			this.id = id;
		}

		public String getId() {
			return id;
		}
	}
	
	private Feature getFeature(String id) {
		for (Feature feature : Feature.values()) {
			if (feature.getId() != null && feature.getId().equals(id)) {
				return feature;
			}
		}
		return Feature.UNRECOGNIZED;
	}
	
	public Exon(Region region, String feature, int exonNumber) {
		this.region = region;
		this.feature = getFeature(feature);
		this.exonNumber = exonNumber;
	}

	public Exon(Region region, String feature, int exonNumber,
			String geneId, String transcId, String geneName, String transcName,
			String biotype) {

		this(region, feature, exonNumber);
		this.geneId = geneId;
		this.transcriptId = transcId;
		this.geneName = geneName;
		this.transcName = transcName;
		this.biotype = biotype;
	}

	public int compareTo(Exon other) {
		
		//All exons here should be from same transcript		
		int transcriptIdComparison = this.transcriptId.compareTo(other.transcriptId);		
		//Sort higher level features (transcript, exon) before lower level features (CDS, start codon),
		//because this is a practical drawing order
		int featureComparison = this.feature.compareTo(other.feature);		
		//Together, transcript, feature and start position should make an unique identifier, so that data isn't lost
		int startComparison = region.start.compareTo(other.region.start);
		
		if (transcriptIdComparison != 0) {
			return transcriptIdComparison;
		} else if (featureComparison != 0){
			return featureComparison;
		} else {
			return startComparison;
		}
	}	
	
	public Object getExonNumber() {
		return exonNumber;
	}

	public Feature getFeature() {
		return feature;
	}

	@Override
	public int hashCode() {
		return transcriptId.hashCode() << 8 + exonNumber << 2 + feature.ordinal();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Exon) {
			return this.compareTo((Exon) o) == 0;
		} else {
			return false;
		}
	}
	
	@Override
	public String toString() {
		
		return region.toString(true) + " \t" + feature + " \t" + transcriptId + " \t" + transcName + " \t" + exonNumber;
	}

	public Region getRegion() {
		return region;
	}

	public Integer getIndex() {
		return exonNumber;
	}

	public void setTranscript(Transcript transc) {
		this.transcript = transc;
	}
	
	public String getGeneId() {
		return geneId;
	}

	public String getTranscriptId() {
		return transcriptId;
	}

	public String getGeneName() {
		return geneName;
	}

	public String getTranscName() {
		return transcName;
	}

	public String getBiotype() {
		return biotype;
	}
}
