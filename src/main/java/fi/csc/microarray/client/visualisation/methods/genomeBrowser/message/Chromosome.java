package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;

public class Chromosome implements Comparable<Chromosome> {
	
	private String chr;
        private Integer intValue;
	
	//Static constants that compareTo method considers bigger and smaller than any real choromosome
	public static final String MAX_VALUE = "ZZZZZZZZ";
	public static final String MIN_VALUE = "0";

	public static final Chromosome MAX = new Chromosome(MAX_VALUE);
	public static final Chromosome MIN = new Chromosome(MIN_VALUE);
	
	public Chromosome(String chr) {
		
		this.chr = chr;
                try {
			intValue = Integer.parseInt(chr);
		} catch (NumberFormatException e) {
		}
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

		if(intValue != null && o.intValue != null) {
			return intValue.compareTo(o.intValue);
		} else if ( intValue != null && o.intValue == null) {
			
			return -1;
			
		} else if ( intValue == null && o.intValue != null) {
			
			return 1;
		} else {
			return chr.compareTo(o.chr);
		}
	}
	
	public String toString() {
		return chr;
	}
}
