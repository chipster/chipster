package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;

public class Chromosome implements Comparable<Chromosome> {
	
	private String chr;
	
	//Static constants that compareTo method considers bigger and smaller than any real choromosome
	public static final String MAX_VALUE = "ZZZZZZZZ";
	public static final String MIN_VALUE = "0";

	public static final Chromosome MAX = new Chromosome(MAX_VALUE);
	public static final Chromosome MIN = new Chromosome(MIN_VALUE);
	
	public Chromosome(String chr) {
		
		this.chr = chr;
	}
	
	public boolean equals(Object o){
		if (o instanceof Chromosome) {
			Chromosome other = (Chromosome) o;
			return chr.equals(other.chr);
		}
		return false;
	}
	
	public int hashCode(){
		return chr.hashCode();
	}

	public int compareTo(Chromosome o){
		
		Integer one = null;
		Integer other = null;
		
		try {		
			one = Integer.parseInt(chr);	
		} catch (NumberFormatException e) {			
		}
		
		try {
			other = Integer.parseInt(o.chr);	
		} catch (NumberFormatException e) {			
		}

		if(one != null && other != null) {
			return one.compareTo(other);
		} else if ( one != null && other == null) {
			
			return -1;
			
		} else if ( one == null && other != null) {
			
			return 1;
		} else {
			return chr.compareTo(o.chr);
		}
	}
	
	public String toString() {
		return chr;
	}
}
