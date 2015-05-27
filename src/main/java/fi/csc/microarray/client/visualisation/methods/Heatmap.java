package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

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
import org.jfree.chart.plot.GradientColorPalette;
import org.jfree.chart.plot.HCPlot;
import org.jfree.chart.plot.HCPlot.Selection;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.hc.HCDataset;
import org.jfree.data.hc.HeatMap;

import fi.csc.microarray.client.selection.SelectionEvent;
import fi.csc.microarray.client.selection.IntegratedSelectionManager;
import fi.csc.microarray.client.visualisation.TableAnnotationProvider;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel.SelectionChangeListener;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.ErrorReportAsException;

/**
 * Classic 2D heatmap of continuous values.
 * 
 * @author Aleksi Kallio
 * 
 */
public class Heatmap extends ChipVisualisation implements PropertyChangeListener, SelectionChangeListener {

	private HCPlot hcPlot;

	protected SelectableChartPanel selectableChartPanel;

	protected Set<Integer> selected = new HashSet<Integer>();

	protected DataBean selectionBean;

	protected JPanel zoomChangerPanel;

	protected JPanel spaceFiller;

	protected JScrollPane scroller;

	protected Dimension preferredSize;

	int rowCount;

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(Heatmap.class);

	@Override
	public JComponent getVisualisation(DataBean data) throws ErrorReportAsException {
		try {

			// Connect selections to correct dataset
			TableAnnotationProvider annotationProvider = new TableAnnotationProvider(data);

			// Create heatmap
			QueryResult heatMapFeature = data.queryFeatures("/column/chip.*");
			LinkedList<String> columns = new LinkedList<String>();
			try (Table heatMapDataIterator = heatMapFeature.asTable()) {

				// Count heatmap rows
				rowCount = 0;
				while (heatMapDataIterator.nextRow()) {
					rowCount++;
				}

				// Count columns that contain expression values
				for (String columnName : heatMapDataIterator.getColumnNames()) {
					columns.add(columnName);
				}
			}
			int columnCount = columns.size();

			HeatMap heatMap = new HeatMap("Heatmap", rowCount, columnCount);

			try (Table heatMapData = data.queryFeatures("/column/*").asTable()) { // fetch all columns to get row names

				int row = -1; // This is increased to 0 in the beginning of the
				// loop

				while (heatMapData.nextRow()) {

					row++;

					String geneName = heatMapData.getStringValue(" ");
					geneName = annotationProvider.getAnnotatedRowname(geneName);
					heatMap.setRowName(row, geneName);

					int i = -1;
					for (String columnName : columns) {

						i++;

						heatMap.update(row, i, heatMapData.getFloatValue(columnName));
					}
				}
			}

			// Set column names
			int i = -1; // increased once before action
			for (String columnName : columns) {

				String sampleName = columnName.substring("chip.".length());
				String realName = data.queryFeatures("/phenodata/linked/describe/" + sampleName).asString();

				i++;

				heatMap.setColumnName(i, realName);
			}

			// Create the chart
			boolean tooltips = true;
			HCDataset dataset = new HCDataset(heatMap, null, null);
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

				this.hcPlot.addChangeListener(new PlotChangeListener() {
					public void plotChanged(PlotChangeEvent event) {
						if (event instanceof ClusteringTreeChangeEvent) {
							// HierarchicalClustering.this.orders.updateVisibleIndexes();
							// HierarchicalClustering.this.updateSelectionsFromApplication(false);
						}
					}
				});

				// Set tooltips
				// if (tooltips) {
				// hcPlot.setToolTipGenerator(new MicroarrayHCToolTipGenerator());
				// }

				// Colors
				double min = getMinValue(dataset.getHeatMap());
				double max = getMaxValue(dataset.getHeatMap());

				GradientColorPalette colors = new GradientColorPalette(new double[] { min, max }, new Color[] { Color.BLUE, Color.BLACK, Color.RED });

				hcPlot.setColoring(colors);

			}

			chart.setTitle((TextTitle) null);

			this.selectableChartPanel = new SelectableChartPanel(chart, this, false);
			this.selectableChartPanel.getChartPanel().addChartMouseListener((HCPlot) chart.getPlot());

			updateSelectionsFromApplication(false);
			application.addClientEventListener(this);

			int blockSize = 10;

			int width = (int) (heatMap.getColumnsCount() * blockSize + hcPlot.getRowTreeSize() + hcPlot.getRowNamesSize() + hcPlot.getLeftMarginSize() + hcPlot.getRightMarginSize());

			// Column tree not visible
			int height = (int) (heatMap.getRowCount() * blockSize + hcPlot.getColumnNamesSize() + hcPlot.getTopMarginSize() + hcPlot.getBottomMarginSize());

			this.preferredSize = new Dimension(width, height);

			this.zoomChangerPanel = new JPanel(new BorderLayout());
			this.spaceFiller = new JPanel();
			((FlowLayout) spaceFiller.getLayout()).setAlignment(FlowLayout.LEFT);
			this.spaceFiller.setBackground(Color.white);
			this.scroller = new JScrollPane(spaceFiller);

			setScaledMode(false);

			return zoomChangerPanel;

		} catch (Exception e) {
			// these are very tricky, mostly caused by bad data
			logger.error(e); // log actual cause
			throw new ErrorReportAsException("Hierarchical clustering cannot be shown.", "The problem is probably caused by unsupported data, such as gene names that have illegal characters in them.", e);
		}
	}

	protected void updateSelectionsFromApplication(boolean dispatchEvent) {
		IntegratedSelectionManager manager = application.getSelectionManager().getSelectionManager(selectionBean);

		selected.clear();
		for (int i : manager.getSelectionAsRows()) {
			selected.add(i);
		}

		showSelection(dispatchEvent);
	}

	public static Double getMinValue(HeatMap heatmap) {
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

	public static Double getMaxValue(HeatMap heatmap) {
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

	protected void showSelection(boolean dispatchEvent) {

		Selection[] detailedSelection = new Selection[rowCount];
		Arrays.fill(detailedSelection, Selection.NO);

//		detailedSelection = calculateRowSelectionDetails();

		hcPlot.showSelection(detailedSelection, true);

//		list.setSelectedRows(selected, this, dispatchEvent, selectionBean);
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

	public void selectionChanged(Rectangle2D.Double selectionRect) {

		if (selectionRect == null) {
			selected.clear();
		} else {

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
//						newSelection.addAll(orders.visibleToBean(entity.getRow()));
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

}



