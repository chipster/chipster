package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel.SelectionChangeListener;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class Volcanoplot extends Scatterplot implements ActionListener, PropertyChangeListener, SelectionChangeListener {

	/**
	 * Analytical expression is -log(p), where p is the first available p.* column.
	 */
	private static final String[] Y_AXIS_COLUMN_HEADERS = new String[] {"p.", "pvalue", "padj", "PValue", "FDR"};
	private static final String[] X_AXIS_COLUMN_HEADERS = new String[] {"FC", "log2FoldChange", "logFC"};

	private float ROUNDING_LIMIT;

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public JPanel createSettingsPanel() {

		JPanel settingsPanel = new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

		xBox = new JComboBox<Variable>();
		yBox = new JComboBox<Variable>();

		useButton = new JButton("Draw");
		useButton.addActionListener(this);
		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 0;
		c.insets.set(10, 10, 10, 10);
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weighty = 0;
		c.weightx = 1.0;
		settingsPanel.add(new JLabel("Fold change "), c);
		c.gridy++;
		settingsPanel.add(xBox, c);
		c.gridy++;
		settingsPanel.add(new JLabel("p-value"), c);
		c.gridy++;
		settingsPanel.add(yBox, c);
		c.gridy++;
		settingsPanel.add(useButton, c);
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1.0;
		settingsPanel.add(new JPanel(), c);

		xBox.addActionListener(this);
		yBox.addActionListener(this);

		return settingsPanel;
	}

	protected void refreshAxisBoxes(DataBean data) {
		if (paramPanel == null) {
			throw new IllegalStateException("must call getParameterPanel first");
		}

		ArrayList<Variable> xCols = new ArrayList<>();
		ArrayList<Variable> yCols = new ArrayList<>();
		
		for (String col : X_AXIS_COLUMN_HEADERS) {
			xCols.addAll(Arrays.asList(VisualisationUtilities.getVariablesFilteredInclusive(data, col, false)));
		}
		for (String col : Y_AXIS_COLUMN_HEADERS) {
			yCols.addAll(Arrays.asList(VisualisationUtilities.getVariablesFilteredInclusive(data, col, false)));
		}
		
		Visualisation.fillComboBox(xBox, xCols.toArray(new Variable[0]));
		Visualisation.fillComboBox(yBox, yCols.toArray(new Variable[0]));
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		this.data = data;

		refreshAxisBoxes(data);

		List<Variable> vars = getFrame().getVariables();

		// If this a redraw from the settings panel, use asked columns
		if (vars != null && vars.size() == 2) {
			xBox.setSelectedItem(vars.get(0));
			yBox.setSelectedItem(vars.get(1));
		}

		xVar = (Variable) xBox.getSelectedItem();
		yVar = (Variable) yBox.getSelectedItem();

		PlotDescription description = new PlotDescription(data.getName(), "fold change (log2)", "-log(p)");

		NumberAxis domainAxis = new NumberAxis(description.xTitle);
		NumberAxis rangeAxis = new NumberAxis(description.yTitle);

		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
		renderer.setLinesVisible(false);
		renderer.setShapesVisible(true);
		renderer.setSeriesPaint(0, Color.green);
		renderer.setSeriesPaint(1, Color.red);
		renderer.setSeriesPaint(2, Color.black);
		renderer.setSeriesPaint(3, Color.lightGray);
		renderer.setShape(new Ellipse2D.Float(-2, -2, 4, 4));

		plot = new XYPlot(new XYSeriesCollection(), domainAxis, rangeAxis, renderer);

		this.updateSelectionsFromApplication(false);
		
		// rounding limit is calculated in updateSelectionsFromApplication
		plot.getRangeAxis().setRange(new Range(0, -Math.log10(ROUNDING_LIMIT)));
		
		JFreeChart chart = new JFreeChart(description.plotTitle, plot);

		chart.removeLegend();

		application.addClientEventListener(this);

		selectableChartPanel = new SelectableChartPanel(chart, this);
		return selectableChartPanel;
	}

	private Iterator<Float> getXValueIterator() throws MicroarrayException {

		return data.queryFeatures(xVar.getExpression()).asFloats().iterator();
	}

	private Iterator<Float> getYValueIterator() throws MicroarrayException {

		return new YValueIterator();
	}

	/**
	 * Class tries to find out the rounding limit of y-values and changes zero values into this limit. Iterated values are also translated
	 * with -log().
	 * 
	 * @author klemela
	 * 
	 */
	private class YValueIterator implements Iterator<Float> {

		private static final float DEFAULT_ROUNDING_LIMIT = 0.001f;
		// "/column/" part of the query comes from the getExpression function
		Iterator<Float> original;

		public YValueIterator() throws MicroarrayException {

			original = data.queryFeatures(yVar.getExpression()).asFloats().iterator();

			// Find smallest non-zero value to find out rounding limit
			float min = Float.MAX_VALUE;

			while (original.hasNext()) {

				float y = original.next();
				if (y < min && y > 0) {
					min = y;
				}
			}

			// Rounding to the nearest 1*10^-n
			// plus one to hide points going into lines because of rounding
			ROUNDING_LIMIT = (float) Math.pow(10, Math.ceil(Math.log10(min)) + 1);
			
			// Sanity check
			if (ROUNDING_LIMIT <= 0 || ROUNDING_LIMIT > DEFAULT_ROUNDING_LIMIT) {
				ROUNDING_LIMIT = DEFAULT_ROUNDING_LIMIT;
			}

			original = data.queryFeatures(yVar.getExpression()).asFloats().iterator();

		}

		public boolean hasNext() {
			return original.hasNext();
		}

		public Float next() {
			float y = original.next();

			if (y < ROUNDING_LIMIT) {
				y = ROUNDING_LIMIT;
			}
			return (float) -Math.log10(y);
		}

		public void remove() {
			original.remove();
		}
	};

	protected void updateXYSerieses() throws MicroarrayException {

		Iterator<Float> xValues = getXValueIterator();
		Iterator<Float> yValues = getYValueIterator();

		XYSeries greenSeries = new XYSeries("");
		XYSeries blackSeries = new XYSeries("");
		XYSeries redSeries = new XYSeries("");
		XYSeries selectedSeries = new XYSeries("");

		int row = 0;
		while (xValues.hasNext() && yValues.hasNext()) {

			float x = xValues.next();
			float y = yValues.next();

			boolean overYThreshold = y >= -Math.log10(0.05);
			boolean overXThreshold = Math.abs(x) >= 1f;

			if (selectedIndexes.contains(row)) {
				selectedSeries.add(new XYDataItem(x, y));
			} else {

				if (overYThreshold && overXThreshold) {
					if (x < 0) {
						greenSeries.add(new XYDataItem(x, y));

					} else {
						redSeries.add(new XYDataItem(x, y));

					}
				} else {
					blackSeries.add(new XYDataItem(x, y));

				}
			}
			row++;
		}

		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(greenSeries);
		dataset.addSeries(redSeries);
		dataset.addSeries(blackSeries);
		dataset.addSeries(selectedSeries);

		plot.setDataset(dataset);

	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return super.canVisualise(bean) && bean.hasTypeTag(MicroarrayModule.TypeTags.SIGNIFICANT_EXPRESSION_FOLD_CHANGES) ;
	}

	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == useButton) {
			List<Variable> vars = new ArrayList<Variable>();
			vars.add((Variable) xBox.getSelectedItem());
			vars.add((Variable) yBox.getSelectedItem());

			application.setVisualisationMethod(new VisualisationMethodChangedEvent(this, MicroarrayModule.VisualisationMethods.VOLCANOPLOT, vars, getFrame().getDatas(), getFrame().getType(), getFrame()));
		}
	}

	public void selectionChanged(Rectangle.Double newSelection) {

		if (newSelection == null) {
			selectedIndexes.clear();
		} else {

			Iterator<Float> xValues;
			Iterator<Float> yValues;

			try {

				xValues = getXValueIterator();
				yValues = getYValueIterator();

				for (int i = 0; xValues.hasNext() && yValues.hasNext(); i++) {

					if (newSelection.contains(new Point.Double(xValues.next(), yValues.next()))) {

						if (selectedIndexes.contains(i)) {
							// Remove from selection if selected twice
							selectedIndexes.remove(i);
						} else {
							selectedIndexes.add(i);
						}
					}
				}
			} catch (MicroarrayException e) {
				application.reportException(e);
			}
		}

		this.list.setSelectedRows(selectedIndexes, this, true, data);

		try {
			updateXYSerieses();
			
		} catch (MicroarrayException e) {
			application.reportException(e);
		}
	}
}
