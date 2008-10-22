package fi.csc.microarray.cluster;

public class ClusterBranchNode extends ClusterNode {
	
	private ClusterNode leftBranch;
	private ClusterNode rightBranch;

	public void setLeftBranch(ClusterNode node) {
		leftBranch = node;
	}

	public void setRightBranch(ClusterNode node) {
		rightBranch = node;
	}

	public ClusterNode getLeftBranch() {
		return leftBranch;
	}

	public ClusterNode getRightBranch() {
		return rightBranch;
	}

	@Override
	public int getLeafCount() {
		return leftBranch.getLeafCount() + rightBranch.getLeafCount();
	}
}
