package fi.csc.microarray.cluster;

public abstract class ClusterNode {
	
	private String length;

	public void setLength(String length) {
		this.length = length;
	}

	public String getLength() {
		return length;
	}

	public abstract int getLeafCount();
}
