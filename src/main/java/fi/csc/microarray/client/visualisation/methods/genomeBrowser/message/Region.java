package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;

public class Region implements Comparable<Region>{
	public Long start;
	public long end;
	
	public Region(long start, long end) {
		this.start = start;
		this.end = end;
	}

	public Region() {
		this(-1, -1);
	}

	public long getLength(){
		return end - start;
	}
	
	public long getMid(){
		return (start + end) / 2; 
	}
	
	
	public String toString(){
		return start + " - " + end;
	}
	
	public Region clone(){
		return new Region(start, end);
	}
	
	public boolean intercepts(Region other){
		return ((this.end > other.start) && this.start < other.end);
	}

	public Region intercept(Region other) {
		return new Region( Math.max(this.start, other.start), Math.min(this.end, other.end));
	}

	public void move(long move) {
		start += move;
		end += move;
	}

	public int compareTo(Region o) {
		return start.compareTo(o.start);
	}
	
	@Override
	public boolean equals(Object o){
		if(o instanceof Region){
			return start.equals(((Region)o).start);
		}
		return false;
	}
	
	@Override
	public int hashCode(){
		return start.hashCode();
	}

	public boolean contains(Long point) {
		return (point >= start && point < end);
	}
}
