package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Chrosome coordinate. 
 */
public class Chromosome implements Comparable<Chromosome> {

	
	private static final String CHROMOSOME_PREFIX = "chr";
	private String chr;
	private Integer intValue;

	public Chromosome(String chr) {
		
//		if (!chr.startsWith(CHROMOSOME_PREFIX)) {
//			chr = CHROMOSOME_PREFIX + chr; 
//		}
		
		this.chr = chr;
		try {
			this.intValue = Integer.parseInt(chr);
		} catch (NumberFormatException e) {
			// was not numeric
			this.intValue = null;
		}
	}

	public Chromosome(Chromosome chr2) {
		this(chr2.chr);
	}

	public boolean equals(Object o) {
		if (o instanceof Chromosome) {
			Chromosome other = (Chromosome) o;
			return chr.equals(other.chr);
		}
		return false;
	}

	public int hashCode() {
		return chr.hashCode();
	}

	public int compareTo(Chromosome o) {

		if (intValue != null && o.intValue != null) {
			return intValue.compareTo(o.intValue);
			
		} else if (intValue != null && o.intValue == null) {
			return -1;

		} else if (intValue == null && o.intValue != null) {
			return 1;
			
		} else {
			return chr.compareTo(o.chr);
		}
	}

	public String toString() {
		return chr;
	}
}
