package fi.csc.microarray.client.visualisation.methods;

import java.awt.event.ActionListener;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.beans.PropertyChangeListener;
import java.util.Iterator;

import javax.swing.JComponent;
import javax.swing.JPanel;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.OverlayLayout;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.databeans.DataBean;

public class Volcanoplot extends Scatterplot implements ActionListener, MouseListener, MouseMotionListener, PropertyChangeListener {

	/**
	 * Analytical expression is -log(p), where p is the first available p.* column.
	 */
	private static final String Y_AXIS_EXPRESSION = "neg(log(/column/p.*))";
	/**
	 * Fold change.
	 */
	private static final String X_AXIS_EXPRESSION = "/column/FC"; 

	public Volcanoplot(VisualisationFrame frame) {
		super(frame);
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		this.data = data;

		// Hide the panel if it exists already ( setting change) to avoid
		// repaints when thes allItems aren't collected yet
		if (chartPanel != null) {
			chartPanel.setVisible(false);
		}

		allItems.clear();

		Iterator<Float> xValues = data.queryFeatures(X_AXIS_EXPRESSION).asFloats().iterator();
		Iterator<Float> yValues = data.queryFeatures(Y_AXIS_EXPRESSION).asFloats().iterator();

		XYSeries redSeries = new XYSeries("", false); // autosort=false, autosort would mess up selection
		XYSeries blackSeries = new XYSeries("", false); // autosort=false, autosort would mess up selection
		int row = 0;
		for (String name : data.queryFeatures("/identifier").asStrings()) {
			float x = xValues.next();
			float y = yValues.next();
			boolean overThresholds = Math.abs(x) >= 2f && y >= 1f;
			int series;
			int index;
			if (overThresholds) {
				series = 1;
				index = blackSeries.getItemCount();
				blackSeries.add(new XYDataItem(x, y));				
			} else {
				series = 0;
				index = redSeries.getItemCount();
				redSeries.add(new XYDataItem(x, y));
			}
			allItems.add(new DataItem2D(null, name, row, index, series));
			row++;
		}

		PlotDescription description = new PlotDescription(data.getName(), "fold change", "-log(p)");

		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(redSeries);
		dataset.addSeries(blackSeries);

		NumberAxis domainAxis = new NumberAxis(description.xTitle);
		domainAxis.setAutoRangeIncludesZero(false);
		NumberAxis rangeAxis = new NumberAxis(description.yTitle);
		rangeAxis.setAutoRangeIncludesZero(false);
		XYPlot plot = new XYPlot(dataset, domainAxis, rangeAxis, null);
		JFreeChart chart = new JFreeChart(description.plotTitle, plot);
		chartPanel = makePanel(chart);
		chartPanel.addMouseListener(this);
		chartPanel.addMouseMotionListener(this);
		chartPanel.setMouseZoomable(false);
		
		// rendered depends on chartPanel for coordinate translations (a bit awkward...), so it cannot be created earlier
		XYItemRenderer renderer = new PositionRecordingRenderer(StandardXYItemRenderer.SHAPES, allItems, selectedItems, chartPanel);
		plot.setRenderer(renderer); 

		overlayPanel = new JPanel(new OverlayLayout());
		transparentPanel = new TransparentPanel();
		overlayPanel.add(transparentPanel);
		overlayPanel.add(chartPanel);

		this.updateSelectionsFromApplication(false);

		application.addPropertyChangeListener(this);

		chartPanel.setVisible(true);

		return overlayPanel;
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		boolean isTabular = VisualisationMethod.SPREADSHEET.getHeadlessVisualiser().canVisualise(bean);
		return isTabular && hasRows(bean) && bean.queryFeatures(Y_AXIS_EXPRESSION).exists() && bean.queryFeatures(X_AXIS_EXPRESSION).exists();
	}
}
