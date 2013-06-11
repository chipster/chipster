package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.SynonymReplace;

/**
 * Utility for converting normalised chromosome names back to original unnormalised format.
 * Unnormaliser initialises itself by taking an unnormalised chromosome name,
 * calling {@link Chromosome#toNormalisedString()} to normalise it, comparing
 * the result to original and storing the prefix and postfix.
 * 
 * @author Aleksi Kallio
 *
 */
public class ChromosomeNameUnnormaliser {
	
	private String prefix;
	private String postfix;
	
	private SynonymReplace synonymReplace = new SynonymReplace();

	public static ChromosomeNameUnnormaliser newIdentityPreversingUnnormaliser() {
		return new ChromosomeNameUnnormaliser("", "");
	}
	
	/**
	 * This constructor accepts a SynonymReplace object for defining format specific synonym versions.   
	 * 
	 * @param prefix
	 * @param postfix
	 * @param synonymReplace
	 */	
	private ChromosomeNameUnnormaliser(String prefix, String postfix) {
		this.prefix = prefix;
		this.postfix = postfix;
	}
	
	public ChromosomeNameUnnormaliser(String unnormalisedName) {
		String normalisedName = Chromosome.normalise(unnormalisedName, false);
		int matchStarts = unnormalisedName.indexOf(normalisedName);		
					
		this.prefix = unnormalisedName.substring(0, matchStarts);
		this.postfix = unnormalisedName.substring(matchStarts + normalisedName.length());
	}
	
	
	/**
	 * Use this constructor when the list of chromosome names is available. The chromosome name
	 * synonyms are search from the list so that unnormalisation later returns those synonyms that 
	 * are used in the file. 
	 * 
	 * @param unnormalisedNames
	 */
	public ChromosomeNameUnnormaliser(List<String> unnormalisedNames) {
		
		this(unnormalisedNames.get(0));
		
		List<String> normalisedNames = new LinkedList<>();
		
		//remove prefix and postfix
		for (String unnormalised : unnormalisedNames) {
			normalisedNames.add(Chromosome.normalise(unnormalised, false));
		}		
				
		for (SynonymReplace.Synonym synonym : Chromosome.getSynonymReplace().getSynonyms()) {
			if (normalisedNames.contains(synonym.getReplaceWith())) {
				//The file contains the default name
			} else if (normalisedNames.contains(synonym.getSearchFor())) {
				synonymReplace.add(new SynonymReplace.Synonym(synonym.getReplaceWith(), synonym.getSearchFor()));
			}
		}
	}

	public String unnormalise(Chromosome normalised) {
		return unnormalise(normalised.toNormalisedString()); 
	}

	public String unnormalise(String normalised) {
		normalised = synonymReplace.apply(normalised);
		return prefix + normalised + postfix;
	}
}
