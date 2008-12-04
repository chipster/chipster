package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

import javax.swing.JComponent;

import org.apache.log4j.Logger;
import org.jfree.chart.BioChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.labels.StandardHCToolTipGenerator;
import org.jfree.chart.plot.GradientColorPalette;
import org.jfree.chart.plot.HCPlot;
import org.jfree.data.hc.DataRange;
import org.jfree.data.hc.DataRangeMismatchException;
import org.jfree.data.hc.HCDataset;
import org.jfree.data.hc.HCTreeNode;
import org.jfree.data.hc.HeatMap;

import fi.csc.microarray.ErrorReportAsException;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.cluster.ClusterBranchNode;
import fi.csc.microarray.cluster.ClusterLeafNode;
import fi.csc.microarray.cluster.ClusterNode;
import fi.csc.microarray.cluster.ClusterParser;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.Table;

public class HierarchicalClustering extends Visualisation {

	public HierarchicalClustering(VisualisationFrame frame) {
		super(frame);
	}

	private static final Logger logger = Logger.getLogger(HierarchicalClustering.class);

	private Map<String, Integer> orderMap = new HashMap<String, Integer>();

	private class MicroarrayHCToolTipGenerator extends StandardHCToolTipGenerator {

		/**
		 * Creates tooltip for Hierarchical clustering visualisation. The
		 * tooltip includes chip and gene names and value.
		 * 
		 * The code is mainly copied from StandardHCToolTipGenerator
		 * 
		 * Note! The rows and columns are in different order than in original
		 * version of this method. In the Viski library the rows and columns are
		 * messed up.
		 */
		@Override
		public String generateToolTip(HeatMap heatmap, DataRange rowRange, DataRange columnRange) {
			int minRow;
			int maxRow;
			int minColumn;
			int maxColumn;
			int rowCounter;
			int columnCounter;
			int blockCount;
			double averageValue;

			try {
				minRow = rowRange.getLeftBound();
				maxRow = rowRange.getRightBound();
				minColumn = columnRange.getLeftBound();
				maxColumn = columnRange.getRightBound();
			} catch (Exception e) {
				return "This block contains no data.";
			}

			if ((minRow == maxRow) && (minColumn == maxColumn)) {
				return "(" + heatmap.getRowName(minRow) + "," + heatmap.getColumnName(minColumn) + ") = " + heatmap.get(minRow, minColumn);
			}

			for (averageValue = 0, blockCount = 0, rowCounter = minRow; rowCounter <= maxRow; rowCounter++) {
				for (columnCounter = minColumn; columnCounter <= maxColumn; columnCounter++, blockCount++) {
					averageValue += heatmap.get(rowCounter, columnCounter);
				}
			}
			averageValue = averageValue / blockCount;

			return "(" + heatmap.getRowName(minRow) + "," + heatmap.getColumnName(minColumn) + ") .. " + "(" + heatmap.getRowName(maxRow) + "," + heatmap.getColumnName(maxColumn) + ") = " + averageValue + " (contains " + (blockCount + 1) + " blocks)";

		}
	}

	/**
	 * getVisualisation method has too modes. If the genes are clustered the
	 * mode is normal and in reversed mode the chips are clustered. Both modes
	 * are implemented in turns of couple rows which makes the structure very
	 * ugly and vulnerable for mistakes. Maybe some day there is time to
	 * restructure this method.
	 * 
	 * @throws MicroarrayException
	 *             is parsing or visualisation generation fails
	 */
	@Override
	public JComponent getVisualisation(DataBean data) throws MicroarrayException {
		try {

			String hcTree = data.queryFeatures("/clusters/hierarchical/tree").asStrings().iterator().next();

			// create heatmap
			QueryResult heatMapFeature = data.queryFeatures("/clusters/hierarchical/heatmap");			
			//TableAnnotationProvider annotationProvider = new TableAnnotationProvider(heatMapFeature.asTable());
			Table heatMapData = heatMapFeature.asTable();
			
			// count rows
			int rowCount = 0;
			while (heatMapData.nextRow()) {
				rowCount++;
			}

			LinkedList<String> columns = new LinkedList<String>();
			for (String columnName : heatMapData.getColumnNames()) {
				if (columnName.startsWith("chip.")) {
					columns.add(columnName);
				} else {
					logger.debug("Column skipped in HC: " + columnName);
				}
			}
			int columnCount = columns.size();

			// HC tree has to be done first to find out the order of rows and number of rows
			HCTreeNode root = null;
			ClusterParser parser = new ClusterParser(hcTree);
			ClusterBranchNode tree = parser.getTree();
			
			
			// check which way we have clustered
			boolean reversed = hcTree.contains("chip.");
			HeatMap heatMap = null;
			if (!reversed) {
				rowCount = tree.getLeafCount(); // heatmap has more genes than tree (sampling done), correct for it
				heatMap = new HeatMap("Heatmap", rowCount, columnCount);
			} else {
				heatMap = new HeatMap("Heatmap", columnCount, rowCount);
			}

			// go through the tree to find its biggest height (depth actually)
			int initialHeight = getTreeHeight(tree);

			orderMap.clear();
			// read the tree and fill the orderMap
			root = readTree(tree, 0, initialHeight);

			// set data to heatmap, reset the iterator
			heatMapData = data.queryFeatures("/clusters/hierarchical/heatmap").asTable();
			int row = -1; // This is increased to 0 in the beginning of the
							// loop

			while (heatMapData.nextRow()) {

				if (!reversed) {
					// Find the row number in heatMap corresponding the name of
					// this row
					// TODO if there is multiple rows with same names we are in
					// trouble
					String key = translate(heatMapData.getStringValue(" "));
					if (orderMap.containsKey(key)) {						
						row = orderMap.get(key);
					} else {
						continue;
					}
					logger.debug("Adding a new row to heatmap, name: " + heatMapData.getStringValue(" ") + "\tto row: " + row);
				} else {
					// reversed row is a column, just use the order from the
					// iteration
					row++;
				}

				String geneName = heatMapData.getStringValue(" ");
				// TODO real gene names cannot be used because the visualisation does not support it
				//geneName = annotationProvider.getAnnotatedRowname(geneName);
				geneName = translate(geneName);

				if (!reversed) {
					heatMap.setRowName(row, translate(heatMapData.getStringValue(" ")));
				} else {
					heatMap.setColumnName(row, translate(heatMapData.getStringValue(" ")));
				}

				int i = -1;
				for (String columnName : columns) {

					if (!reversed) {
						// column index, just use the order from the iteration
						i++;
					} else {
						i = orderMap.get(columnName);
						logger.debug("Adding a new row to heatmap (reversed), name: " + columnName + "\tto row: " + i);
					}

					if (!reversed) {
						heatMap.update(row, i, heatMapData.getFloatValue(columnName));
					} else {
						heatMap.update(i, row, heatMapData.getFloatValue(columnName));
					}
				}
			}

			// Set column names (row names if reversed)
			int i = -1; // increased once before action
			for (String columnName : columns) {

				String sampleName = columnName.substring("chip.".length());
				String realName = data.queryFeatures("/phenodata/linked/sample_to_name/" + sampleName).asString();

				if (!reversed) {
					// column index, just use the order from the iteration
					i++;
				} else {
					i = orderMap.get(columnName);
					logger.debug("Adding a new row to heatmap (reversed), name: " + columnName + "\tto row: " + i);
				}

				if (!reversed) {
					heatMap.setColumnName(i, realName);
				} else {
					heatMap.setRowName(i, realName);
				}
			}

			HCDataset dataset = new HCDataset(heatMap, root, null);

			// create the chart...
			boolean tooltips = true;
			JFreeChart chart = BioChartFactory.createHCChart("Hierarchical Clustering", // chart
																						// title
			dataset, // data
			tooltips, // tooltips?
			false // URLs?
			);

			// set special tooltips to hcChart
			if (chart.getPlot() instanceof HCPlot) {
				HCPlot hcPlot = (HCPlot) chart.getPlot();

				// Add more space to column and row names
				hcPlot.setColumnNamesSize(0.2);
				hcPlot.setRowNamesSize(0.3);

				// Space for the trees
				hcPlot.setRowTreeSize(0.1);
				hcPlot.setColumnTreeSize(0.1);

				// Set tooltips
				if (tooltips) {
					hcPlot.setToolTipGenerator(new MicroarrayHCToolTipGenerator());
				}

				// Colors
				double min = getMinValue(dataset.getHeatMap());
				double max = getMaxValue(dataset.getHeatMap());
				double mid = ((min + max) / 2);

				GradientColorPalette colors = new GradientColorPalette();
				colors.addKeyColor(min, Color.GREEN);
				colors.addKeyColor(mid, Color.BLACK);
				colors.addKeyColor(max, Color.RED);

				hcPlot.setColoring(colors);

			}

			ChartPanel chartPanel = makePanel(chart);
			chartPanel.addChartMouseListener((HCPlot) chart.getPlot());
			return chartPanel;

		} catch (Exception e) {
			// these are very tricky, mostly caused by bad data
			logger.error(e); // log actual cause
			throw new ErrorReportAsException("Hierarchical clustering cannot be shown.", "The problem is probably caused by unsupported data, such as gene names that have illegal characters in them.", e);
		}

	}

	/**
	 * Translates name to parenthesis tree formate. This translation is required
	 * because R script does it.
	 */
	private String translate(String gene) {
		return gene.replace("(", "").replace(")", "-");
	}

	private int getTreeHeight(ClusterNode tree) {
		if (tree instanceof ClusterLeafNode) {
			return 0;
		} else {
			int left = getTreeHeight(((ClusterBranchNode) tree).getLeftBranch()) + 1;
			int right = getTreeHeight(((ClusterBranchNode) tree).getRightBranch()) + 1;

			return left > right ? left : right;
		}
	}

	private HCTreeNode readTree(ClusterNode tree, int index, int height) throws DataRangeMismatchException {
		if (tree instanceof ClusterLeafNode) {
			String gene = ((ClusterLeafNode) tree).getGene();
			orderMap.put(gene, index);
			logger.debug("LeafNode: " + ((ClusterLeafNode) tree).getGene() + "in index: " + index);
			return new HCTreeNode(0, index); // height is zero

		} else {
			HCTreeNode node = new HCTreeNode(height);
			HCTreeNode leftChild = readTree(((ClusterBranchNode) tree).getLeftBranch(), index, height - 1);
			node.setLeftChild(leftChild);
			HCTreeNode rightChild = readTree(((ClusterBranchNode) tree).getRightBranch(), node.getDataRange().getRightBound() + 1, height - 1);
			node.setRightChild(rightChild);

			return node;
		}
	}

	private static Double getMinValue(HeatMap heatmap) {
		Double min = null;
		for (int row = 0; row < heatmap.getRowCount(); row++) {
			for (int column = 0; column < heatmap.getColumnsCount(); column++) {
				Double value = heatmap.get(row, column);
				if (min == null || value < min) {
					min = value;
				}
			}
		}
		return min;
	}

	private static Double getMaxValue(HeatMap heatmap) {
		Double max = null;
		for (int row = 0; row < heatmap.getRowCount(); row++) {
			for (int column = 0; column < heatmap.getColumnsCount(); column++) {
				Double value = heatmap.get(row, column);
				if (max == null || value > max) {
					max = value;
				}
			}
		}
		return max;
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return bean.isContentTypeCompatitible("application/x-treeview") && bean.queryFeatures("/clusters/hierarchical/tree").exists() && bean.queryFeatures("/clusters/hierarchical/heatmap").exists();
	}
}
