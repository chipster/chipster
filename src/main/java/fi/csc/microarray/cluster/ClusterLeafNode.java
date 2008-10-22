package fi.csc.microarray.cluster;

public class ClusterLeafNode extends ClusterNode {

	private String gene;

	public void setGene(String gene) {
		this.gene = gene;
	}

	public String getGene() {
		return gene;
	}

	@Override
	public int getLeafCount() {
		return 1;
	}
}
