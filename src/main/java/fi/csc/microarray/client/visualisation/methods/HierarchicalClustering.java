package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;

import org.apache.log4j.Logger;
import org.jfree.chart.BioChartFactory;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.entity.EntityCollection;
import org.jfree.chart.entity.HCTreeNodeEntity;
import org.jfree.chart.entity.HeatMapBlockEntity;
import org.jfree.chart.event.ClusteringTreeChangeEvent;
import org.jfree.chart.event.PlotChangeEvent;
import org.jfree.chart.event.PlotChangeListener;
import org.jfree.chart.labels.StandardHCToolTipGenerator;
import org.jfree.chart.plot.GradientColorPalette;
import org.jfree.chart.plot.HCPlot;
import org.jfree.chart.plot.HCPlot.Selection;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.hc.DataRange;
import org.jfree.data.hc.DataRangeMismatchException;
import org.jfree.data.hc.HCDataset;
import org.jfree.data.hc.HCTreeNode;
import org.jfree.data.hc.HeatMap;

import fi.csc.microarray.client.selection.SelectionEvent;
import fi.csc.microarray.client.selection.IntegratedSelectionManager;
import fi.csc.microarray.client.visualisation.SelectionList;
import fi.csc.microarray.client.visualisation.TableAnnotationProvider;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel.SelectionChangeListener;
import fi.csc.microarray.client.visualisation.methods.hc.OrderSuperviser;
import fi.csc.microarray.cluster.ClusterBranchNode;
import fi.csc.microarray.cluster.ClusterLeafNode;
import fi.csc.microarray.cluster.ClusterNode;
import fi.csc.microarray.cluster.ClusterParser;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.Traversal;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.ErrorReportAsException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class HierarchicalClustering extends Visualisation implements PropertyChangeListener, SelectionChangeListener {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	protected SelectableChartPanel selectableChartPanel;

	private static final Logger logger = Logger.getLogger(HierarchicalClustering.class);

	OrderSuperviser orders;

	private JPanel paramPanel;
	private SelectionList list;

	// Selected indexes in the order of parent data bean
	protected Set<Integer> selected = new HashSet<Integer>();

	protected DataBean selectionBean;

	private HCPlot hcPlot;

	private boolean reversed;

	protected JPanel zoomChangerPanel;

	protected JPanel spaceFiller;

	protected JScrollPane scroller;

	protected Dimension preferredSize;

	private JCheckBox zoomCheckBox;

	private class MicroarrayHCToolTipGenerator extends StandardHCToolTipGenerator {

		/**
		 * Creates tooltip for Hierarchical clustering visualisation. The tooltip includes chip and gene names and value.
		 * 
		 * The code is mainly copied from StandardHCToolTipGenerator
		 * 
		 * Note! The rows and columns are in different order than in original version of this method. In the Viski library the rows and
		 * columns are messed up.
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

	@Override
	public JPanel getParameterPanel() {
		if (paramPanel == null) {
			paramPanel = new JPanel();
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);
			paramPanel.setLayout(new BorderLayout());

			JPanel settings = this.createSettingsPanel();
			list = new SelectionList();

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Settings", settings);
			tabPane.addTab("Selected", list);

			paramPanel.add(tabPane, BorderLayout.CENTER);
		}
		return paramPanel;
	}

	public JPanel createSettingsPanel() {

		JPanel settingsPanel = new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

		zoomCheckBox = new JCheckBox("Fit to screen", false);

		zoomCheckBox.addActionListener(new ActionListener() {

			public void actionPerformed(ActionEvent e) {
				setScaledMode(zoomCheckBox.isSelected());
			}
		});

		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 0;
		c.insets.set(10, 10, 0, 10);
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weighty = 0;
		c.weightx = 1.0;
		settingsPanel.add(zoomCheckBox, c);
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1.0;
		settingsPanel.add(new JPanel(), c);

		return settingsPanel;
	}

	/**
	 * getVisualisation method has too modes. If the genes are clustered the mode is normal and in reversed mode the chips are clustered.
	 * Both modes are implemented in turns of couple rows which makes the structure very ugly and vulnerable for mistakes. Maybe some day
	 * there is time to restructure this method.
	 * 
	 * @throws MicroarrayException
	 *             is parsing or visualisation generation fails
	 */
	@Override
	public JComponent getVisualisation(DataBean data) throws MicroarrayException {

		try {

			// First find a dataset to which user's selections are connected to.
			// There is no point to connect selections to cluster tree data because it's not possible to visualise it
			// with other visualisation methods. The parent that contains the expression values is used instead.
			List<DataBean> selectionBeans = data.traverseLinks(new Link[] { Link.DERIVATION }, Traversal.DIRECT);

			// First one is the correct one
			if (selectionBeans.size() > 1) {
				selectionBean = selectionBeans.get(1);
			}

			if (selectionBean == null) {
				throw new ErrorReportAsException("Source dataset not found", "Hierarchical clustering " + "needs its source dataset.", " Select both HC and its source dataset by keeping \n" + "Ctrl key pressed and right click with mouse over one of them to create \n" + "derivation link from the original dataset to \n" + "clustered one.");
			}

			// Connect selections to correct dataset
			TableAnnotationProvider annotationProvider = new TableAnnotationProvider(selectionBean);

			// Create heatmap
			QueryResult heatMapFeature = data.queryFeatures("/clusters/hierarchical/heatmap");
			LinkedList<String> columns = new LinkedList<String>();
			int rowCount = 0;
			try (Table heatMapDataIterator = heatMapFeature.asTable()) {

				// Count heatmap rows
				while (heatMapDataIterator.nextRow()) {
					rowCount++;
				}

				// Count columns that contain expression values
				for (String columnName : heatMapDataIterator.getColumnNames()) {
					if (columnName.startsWith("chip.")) {
						columns.add(columnName);
					} else {
						logger.debug("Column skipped in HC: " + columnName);
					}
				}
			}
			int columnCount = columns.size();

			// Parse HC tree and check which way we have clustered
			String hcTree = data.queryFeatures("/clusters/hierarchical/tree").asStrings().iterator().next();
			ClusterBranchNode tree = new ClusterParser(hcTree).getTree();
			this.reversed = hcTree.contains("chip.");

			// Adjust gene count for sampling
			HeatMap heatMap = null;
			int dataCount;
			if (!reversed) {
				rowCount = tree.getLeafCount(); // heatmap has more genes than tree (sampling done), correct for it
				heatMap = new HeatMap("Heatmap", rowCount, columnCount);
				dataCount = rowCount;
			} else {
				heatMap = new HeatMap("Heatmap", columnCount, rowCount);
				dataCount = columnCount;
			}

			// Go through the tree to find its biggest height
			int initialHeight = getTreeHeight(tree);

			// Read the tree and fill the treeToId map
			List<String> treeToId = new ArrayList<String>();
			treeToId.addAll(Collections.nCopies(dataCount, (String) null));
			HCTreeNode root = readTree(tree, 0, initialHeight, treeToId);

			orders = new OrderSuperviser();
			orders.setTreeToId(treeToId);

			List<Integer> treeToBean = new ArrayList<Integer>();
			try (Table heatMapData = data.queryFeatures("/clusters/hierarchical/heatmap").asTable()) {

				treeToBean.addAll(Collections.nCopies(rowCount, -1));

				int row = -1; // This is increased to 0 in the beginning of the
				// loop
				int originalRow = 0;

				while (heatMapData.nextRow()) {

					if (!reversed) {
						// Find the row number in heatMap corresponding the name of
						// this row

						String key = translate(heatMapData.getStringValue(" "));
						if (orders.idToTree(key) != -1) { // if the id is found
							row = orders.idToTree(key);
							treeToBean.set(row, originalRow);
							originalRow++;
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
					geneName = annotationProvider.getAnnotatedRowname(geneName);

					if (!reversed) {
						heatMap.setRowName(row, geneName);
					} else {
						heatMap.setColumnName(row, geneName);
					}

					int i = -1;
					for (String columnName : columns) {

						if (!reversed) {
							// column index, just use the order from the iteration
							i++;
						} else {
							i = orders.idToTree(columnName);
							logger.debug("Adding a new row to heatmap (reversed), name: " + columnName + "\tto row: " + i);
						}

						if (!reversed) {
							heatMap.update(row, i, heatMapData.getFloatValue(columnName));
						} else {
							heatMap.update(i, row, heatMapData.getFloatValue(columnName));
						}
					}
				}
			}

			orders.setTreeToBean(treeToBean);

			// Set column names (row names if reversed)
			int i = -1; // increased once before action
			for (String columnName : columns) {

				String sampleName = columnName.substring("chip.".length());
				String realName = data.queryFeatures("/phenodata/linked/describe/" + sampleName).asString();

				if (!reversed) {
					// column index, just use the order from the iteration
					i++;
				} else {
					i = orders.idToTree(columnName);
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

				this.hcPlot = hcPlot;
				orders.setPlot(hcPlot);

				this.hcPlot.addChangeListener(new PlotChangeListener() {
					public void plotChanged(PlotChangeEvent event) {
						if (event instanceof ClusteringTreeChangeEvent) {
							HierarchicalClustering.this.orders.updateVisibleIndexes();
							HierarchicalClustering.this.updateSelectionsFromApplication(false);
						}
					}
				});

				// Set tooltips
				if (tooltips) {
					hcPlot.setToolTipGenerator(new MicroarrayHCToolTipGenerator());
				}

				// Colors
				double min = Heatmap.getMinValue(dataset.getHeatMap());
				double max = Heatmap.getMaxValue(dataset.getHeatMap());

				GradientColorPalette colors = new GradientColorPalette(new double[] { min, max }, new Color[] { Color.BLUE, Color.BLACK, Color.RED });

				hcPlot.setColoring(colors);

			}

			chart.setTitle((TextTitle) null);

			selectableChartPanel = new SelectableChartPanel(chart, this, false);
			selectableChartPanel.getChartPanel().addChartMouseListener((HCPlot) chart.getPlot());

			updateSelectionsFromApplication(false);
			application.addClientEventListener(this);

			int blockSize = 10;

			int width = (int) (heatMap.getColumnsCount() * blockSize + hcPlot.getRowTreeSize() + hcPlot.getRowNamesSize() + hcPlot.getLeftMarginSize() + hcPlot.getRightMarginSize());

			// Column tree not visible
			int height = (int) (heatMap.getRowCount() * blockSize + hcPlot.getColumnNamesSize() + hcPlot.getTopMarginSize() + hcPlot.getBottomMarginSize());

			preferredSize = new Dimension(width, height);

			zoomChangerPanel = new JPanel(new BorderLayout());
			spaceFiller = new JPanel();
			((FlowLayout) spaceFiller.getLayout()).setAlignment(FlowLayout.LEFT);
			spaceFiller.setBackground(Color.white);
			scroller = new JScrollPane(spaceFiller);

			setScaledMode(false);

			return zoomChangerPanel;

		} catch (Exception e) {
			// these are very tricky, mostly caused by bad data
			logger.error(e); // log actual cause
			throw new ErrorReportAsException("Hierarchical clustering cannot be shown.", "The problem is probably caused by unsupported data, such as gene names that have illegal characters in them.", e);
		}

	}

	public void setScaledMode(boolean scaled) {

		/*
		 * Ugly way to change zoom level by changing containing panel layout and scroller existence, but JFreeChart scaling is little bit
		 * problematic in this kind of usage.
		 */
		if (scaled) {
			spaceFiller.remove(selectableChartPanel);
			zoomChangerPanel.remove(scroller);
			zoomChangerPanel.add(selectableChartPanel, BorderLayout.CENTER);
			selectableChartPanel.setPreferredSize(null);
		} else {
			spaceFiller.add(selectableChartPanel);
			zoomChangerPanel.remove(selectableChartPanel);
			zoomChangerPanel.add(scroller, BorderLayout.CENTER);
			selectableChartPanel.setPreferredSize(preferredSize);
		}

		zoomChangerPanel.validate();
		zoomChangerPanel.repaint();
	}

	/**
	 * Translates name to parenthesis tree format. This translation is required because R script does it.
	 */
	private String translate(String gene) {

		while (gene.startsWith(" ") || gene.startsWith("(")) {
			gene = gene.substring(1);
		}
		while (gene.endsWith(" ") || gene.endsWith(")")) {
			gene = gene.substring(0, gene.length() - 1);
		}

		gene = gene.replace("(", "-").replace(")", "-").replace(" ", "_");

		while (gene.contains("--")) {
			gene = gene.replace("--", "-");
		}
		return gene;
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

	private HCTreeNode readTree(ClusterNode tree, int index, int height, List<String> treeToId) throws DataRangeMismatchException {

		if (tree instanceof ClusterLeafNode) {
			String gene = ((ClusterLeafNode) tree).getGene();
			treeToId.set(index, gene);
			logger.debug("LeafNode: " + ((ClusterLeafNode) tree).getGene() + "in index: " + index);
			return new HCTreeNode(0, index); // height is zero

		} else {
			HCTreeNode node = new HCTreeNode(height);
			HCTreeNode leftChild = readTree(((ClusterBranchNode) tree).getLeftBranch(), index, height - 1, treeToId);

			node.setLeftChild(leftChild);

			HCTreeNode rightChild = readTree(((ClusterBranchNode) tree).getRightBranch(), node.getDataRange().getRightBound() + 1, height - 1, treeToId);

			node.setRightChild(rightChild);

			return node;
		}
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		DataBean parentBean = MicroarrayModule.getProperSource(bean); 
		return bean.isContentTypeCompatitible("application/x-treeview") && parentBean != null && parentBean.hasTypeTag(MicroarrayModule.TypeTags.NORMALISED_EXPRESSION_VALUES);
	}

	public void selectionChanged(Rectangle2D.Double selectionRect) {

		if (selectionRect == null) {
			selected.clear();
		} else {

			orders.updateVisibleIndexes();

			ChartRenderingInfo info = selectableChartPanel.getChartPanel().getChartRenderingInfo();

			EntityCollection entities = info.getEntityCollection();

			Set<Integer> newSelection = new HashSet<Integer>();

			for (Object obj : entities.getEntities()) {

				// Don't clear the selection if tree was clicked
				if (obj instanceof HCTreeNodeEntity) {
					HCTreeNodeEntity entity = (HCTreeNodeEntity) obj;
					if (entity.getArea().intersects(selectionRect)) {
						return;
					}
				}

				if (obj instanceof HeatMapBlockEntity) {
					HeatMapBlockEntity entity = (HeatMapBlockEntity) obj;

					if (entity.getArea().intersects(selectionRect)) {

						if (!reversed) {
							newSelection.addAll(orders.visibleToBean(entity.getRow()));
						} else {
							newSelection.add(entity.getColumn());
						}
					}
				}
			}

			// New selections can't be put directly to selectedIndexes, because every other occurrence
			// of block (several in one line) inside selection rectangle would undo selection

			for (Integer row : newSelection) {
				if (selected.contains(row)) {
					selected.remove(row);
				} else {
					selected.add(row);
				}
			}

			showSelection(true);
		}
	}

	public void propertyChange(PropertyChangeEvent evt) {
		if (evt instanceof SelectionEvent && evt.getSource() != this && ((SelectionEvent) evt).getData() == selectionBean) {
			updateSelectionsFromApplication(false);
		}
	}

	protected void updateSelectionsFromApplication(boolean dispatchEvent) {
		IntegratedSelectionManager manager = application.getSelectionManager().getSelectionManager(selectionBean);

		orders.updateVisibleIndexes();

		selected.clear();
		for (int i : manager.getSelectionAsRows()) {
			selected.add(i);
		}

		showSelection(dispatchEvent);
	}

	protected void showSelection(boolean dispatchEvent) {

		Selection[] detailedSelection;

		if (!reversed) {
			detailedSelection = calculateRowSelectionDetails();
		} else {
			detailedSelection = calculateColumnSelectionDetails();
		}

		hcPlot.showSelection(detailedSelection, !reversed);

		list.setSelectedRows(selected, this, dispatchEvent, selectionBean);
	}

	private Selection[] calculateRowSelectionDetails() {
		// Number of rows represented by each visible row
		int[] closedRows = orders.getCountOfVisibleReferences();
		// Number of selected in each visible row
		int[] selectedRows = new int[closedRows.length];
		// Is each visible row fully, partially or not at all selected
		Selection[] detailedSelection = new Selection[closedRows.length];

		for (int selectedRow : selected) {
			selectedRows[orders.beanToVisible(selectedRow)]++;
		}

		for (int i = 0; i < selectedRows.length; i++) {
			if (selectedRows[i] == 0) {
				detailedSelection[i] = Selection.NO;

			} else if (selectedRows[i] == closedRows[i]) {
				detailedSelection[i] = Selection.YES;

			} else if (selectedRows[i] < closedRows[i]) {
				detailedSelection[i] = Selection.PARTIAL;
			}
		}

		return detailedSelection;
	}

	private Selection[] calculateColumnSelectionDetails() {
		// Columns are already in right order, so just convert the selections to right format
		Selection[] detailedSelection = new Selection[hcPlot.getDataset().getHeatMap().getColumnsCount()];

		for (int i = 0; i < detailedSelection.length; i++) {
			if (selected.contains(i)) {
				detailedSelection[i] = Selection.YES;
			} else {
				detailedSelection[i] = Selection.NO;
			}
		}

		return detailedSelection;
	}
}
