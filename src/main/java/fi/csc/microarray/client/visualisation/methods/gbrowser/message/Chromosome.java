package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * <p>Chrosome part of a genomic coordinate.</p> 
 * 
 * <p>The class supports different naming 
 * conventions: X, chrX, chrX.fa, X.fa. Internally names
 * are normalised and normalised name is used in all comparisons,
 * so naming conventions can be mixed.</p>
 * 
 * @author Petri Klemelä, Aleksi Kallio
 */
public class Chromosome implements Comparable<Chromosome> {

	public static final String CHROMOSOME_PREFIX = "chr";
	private String chr;
	private Integer intValue;
	private String chrNormalised;
	
	//Default
	private static final SynonymReplace synonymReplace = new SynonymReplace( new SynonymReplace.Synonym[] { new SynonymReplace.Synonym("M", "MT")});
	
	//This is useful with files that don't follow Ensembl convention
	public static final SynonymReplace ucscSynonymReplace = new SynonymReplace( new SynonymReplace.Synonym[] { new SynonymReplace.Synonym("MT", "M")});

	public Chromosome(String chr) {
		
		// store original value
		this.chr = chr;
		this.chrNormalised = normalise(chr, true);
		
		try {
			this.intValue = Integer.parseInt(chrNormalised);
		} catch (NumberFormatException e) {
			// was not numeric
			this.intValue = null;
		}
	}

	public static String normalise(String original, boolean enableSynonymReplace) {
		// If contains any postfix, remove it
		if (original.indexOf(".") != -1) {
			original = original.substring(0, original.indexOf("."));
		}
		
		// Remove known prefix, if exists
		String normalised = original.replace(CHROMOSOME_PREFIX, "");
		
		if (enableSynonymReplace) {
			normalised = synonymReplace.apply(normalised);
		}
		
		return normalised;
	}

	public Chromosome(Chromosome chromosome) {
		this(chromosome.chr);
	}

	public boolean equals(Object o) {
		if (o instanceof Chromosome) {
			Chromosome other = (Chromosome) o;
			return chrNormalised.equals(other.chrNormalised);
		}
		return false;
	}

	public int hashCode() {
		return chrNormalised.hashCode();
	}

	/**
     * <p>Compares this chromosome with the specified chromosome for order.  Returns a
     * negative integer, zero, or a positive integer as this object is less
     * than, equal to, or greater than the specified object. Chromosomes
     * that are number are compated numerically. Chromosome that
     * are non-numerical are compared lexically (in their normalised form).
     * Numerical names are always considered smaller than non-numerical.</p>
     * 
     * <p>In the case of human genome, the resulting order is: 1, 2, 3, ..., 21, 22, X, Y.
     * It is the standard order commonly used, which was also the design
     * guideline for this function.</p>
     */
	public int compareTo(Chromosome o) {

		if (intValue != null && o.intValue != null) {
			return intValue.compareTo(o.intValue);
			
		} else if (intValue != null && o.intValue == null) {
			return -1;

		} else if (intValue == null && o.intValue != null) {
			return 1;
			
		} else {
			return chrNormalised.compareTo(o.chrNormalised);
		}
	}

	/**
	 * Returns the original chromosome name.
	 */
	public String getOriginalName() {
		return chr;
	}	
	
	/**
	 * Returns the normalised chromosome name.
	 */
	public String toString() {
		return chrNormalised;
	}

	/**
	 * Returns the normalised chromosome name.
	 */
	public String toNormalisedString() {
		return chrNormalised;
	}
	
	public Chromosome clone() {
		return new Chromosome(this.chr);
	}

	public static SynonymReplace getSynonymReplace() {
		return synonymReplace;
	}
}
