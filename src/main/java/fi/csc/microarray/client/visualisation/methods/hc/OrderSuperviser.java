package fi.csc.microarray.client.visualisation.methods.hc;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.jfree.chart.plot.HCPlot;
import org.jfree.chart.plot.HCTreeNodeInfo;

public class OrderSuperviser {

	HCPlot hcPlot;

	public void setPlot(HCPlot hcPlot) {
		this.hcPlot = hcPlot;
	}

	// M A P S ////////////////////////////////////////////////////////////////

	private List<String> treeToId = new ArrayList<String>();

	// List indexes referring to full tree order, Integer objects referring to original indexes
	private List<Integer> treeToBean = new ArrayList<Integer>();

	// List indexes referring to full tree order, Integer objects referring to visible indexes
	private List<Integer> treeToVisible = new ArrayList<Integer>();

	// S E T T E R S //////////////////////////////////////////////////////////

	public void setTreeToId(List<String> treeToId) {
		this.treeToId = treeToId;
	}

	public void setTreeToBean(List<Integer> treeToBean) {
		this.treeToBean = treeToBean;
	}

	public void setTreeToVisible(List<Integer> treeToVisible) {
		this.treeToVisible = treeToVisible;
	}

	// C O N V E R T E R S ////////////////////////////////////////////////////

	public int idToTree(String id) {
		return treeToId.indexOf(id);
	}

	public int treeToVisible(int tree) {
		return treeToVisible.get(tree);
	}

	/**
	 * The result can be more than one row, because closed rows represent several individual rows
	 * 
	 * @param visibleToFind
	 * @return
	 */
	public List<Integer> visibleToTree(int visibleToFind) {
		List<Integer> treeIndexes = new ArrayList<Integer>();

		int tree = 0;
		for (int visible : treeToVisible) {

			if (visible == visibleToFind) {
				treeIndexes.add(tree);
			}

			tree++;
		}
		return treeIndexes;
	}

	public int treeToBean(int tree) {
		return treeToBean.get(tree);
	}

	public int beanToTree(int bean) {
		return treeToBean.indexOf(bean);
	}

	public List<Integer> visibleToBean(Integer index) {
		List<Integer> beanIndexes = new ArrayList<Integer>();

		for (int tree : visibleToTree(index)) {
			beanIndexes.add(treeToBean(tree));
		}
		return beanIndexes;
	}

	public Integer beanToVisible(Integer bean) {
		return treeToVisible(beanToTree(bean));
	}

	// K E E P D A T A U P T O D A T E ////////////////////////////////////

	public void updateVisibleIndexes() {
		List<Integer> treeToVisible = new ArrayList<Integer>();

		boolean openRoot = hcPlot.getRowClusteringInfo().getRootNode().isNodeOpen();

		updateVisibleIndexesSubTree(hcPlot.getRowClusteringInfo().getRootNode(), treeToVisible, openRoot, false);

		this.setTreeToVisible(treeToVisible);
	}

	/**
	 * Method updateVisibleIndexes calls this to go trough the tree recursively. There should be no need to call this directly.
	 * 
	 * @param node
	 * @param list
	 * @param openSubTree
	 * @param firstLeafOfClosed
	 * @return
	 */
	private boolean updateVisibleIndexesSubTree(HCTreeNodeInfo node, List<Integer> list, boolean openSubTree, boolean firstLeafOfClosed) {

		/*
		 * It's not very easy job to calculate the relations between indexes of partly closed and open tree. Here is one recursive
		 * implementation which seems to work right with tested trees.
		 * 
		 * OpenSubTree is used to carry information if we are in closed part of tree or not, because one closed branch hides also its sub
		 * branches, whereas their internal state still shows the situation when the tree is open.
		 * 
		 * When inside this kind of closed branch, all the leafs get some visible index as the first one of that closed branch. However, the
		 * first has to be still increased by one from the previous row, because it's still a new row in the visualisation. This situation
		 * is recognised with variable firstLeafOfClosed.
		 */

		if ((node.getLeftChild() == null) && (node.getRightChild() == null)) {
			if (list.size() == 0) {
				// First one is always at index zero (and one object is needed in list to calculate
				// next ones)
				list.add(0);

				if (!openSubTree) {
					// Without this the indexes go wrong if first rows is inside closed
					// branch
					firstLeafOfClosed = false;
				}
			} else {
				if (openSubTree || firstLeafOfClosed) {
					firstLeafOfClosed = false;

					// Visible index of open row is index of previous row increased by one
					list.add(list.get(list.size() - 1) + 1);

				} else {
					// Visible index of closed row is same as previous
					list.add(list.get(list.size() - 1));

				}
			}
		} else if (openSubTree && !node.isNodeOpen()) {
			// Beginning of closed branch, openSubTree gets the value of false
			updateVisibleIndexesSubTree(node, list, false, true);

		} else {

			// Just open the left and right children with the same state information
			if (node.getLeftChild() != null)
				firstLeafOfClosed = updateVisibleIndexesSubTree(node.getLeftChild(), list, openSubTree, firstLeafOfClosed);

			if (node.getRightChild() != null)
				firstLeafOfClosed = updateVisibleIndexesSubTree(node.getRightChild(), list, openSubTree, firstLeafOfClosed);
		}

		return firstLeafOfClosed;
	}

	// U T I L S /////////////////////////////////

	public int[] getCountOfVisibleReferences() {

		List<Integer> countOfVisibleReferences = new ArrayList<Integer>();
		countOfVisibleReferences.addAll(Collections.nCopies(treeToVisible.size(), 0));

		int maxVisible = 0;

		for (int visible : treeToVisible) {
			countOfVisibleReferences.set(visible, countOfVisibleReferences.get(visible) + 1);

			if (visible > maxVisible) {
				maxVisible = visible;
			}
		}

		int[] countTable = new int[maxVisible + 1];

		for (int i = 0; i <= maxVisible; i++) {
			countTable[i] = countOfVisibleReferences.get(i);
		}

		return countTable;
	}
}
