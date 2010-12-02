package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;

public class ChromosomeNameUnnormaliser {
	
	private String prefix;
	private String postfix;

	public static ChromosomeNameUnnormaliser newIdentityPreversingUnnormaliser() {
		return new ChromosomeNameUnnormaliser("", "");
	}
	
	private ChromosomeNameUnnormaliser(String prefix, String postfix) {
		this.prefix = prefix;
		this.postfix = postfix;
	}
	
	public ChromosomeNameUnnormaliser(String unnormalisedName) {
		String normalisedName = new Chromosome(unnormalisedName).toNormalisedString();
		int matchStarts = unnormalisedName.indexOf(normalisedName);
		
		this.prefix = unnormalisedName.substring(0, matchStarts);
		this.postfix = unnormalisedName.substring(matchStarts + normalisedName.length());
	}
	
	public String unnormalise(Chromosome normalised) {
		return prefix + normalised.toNormalisedString() + postfix; 
	}

}
