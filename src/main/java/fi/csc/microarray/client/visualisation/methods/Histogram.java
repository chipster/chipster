package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSpinner;
import javax.swing.SpinnerNumberModel;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.labels.ItemLabelAnchor;
import org.jfree.chart.labels.ItemLabelPosition;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.BarRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;
import org.jfree.ui.TextAnchor;

import fi.csc.microarray.client.visualisation.ChipVisualisation;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.FloatArrayList;

public class Histogram extends ChipVisualisation {

	public Histogram(VisualisationFrame frame) {
		super(frame);
	}

	private JPanel paramPanel;
	private JPanel visualisationPanel = new JPanel(new BorderLayout());
	private JSpinner barCountSpinner = new JSpinner(new SpinnerNumberModel(20, 1, Integer.MAX_VALUE, 1));
	private JButton updateButton;
	private JComboBox chipBox;
	private DataBean data;

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		this.data = data;
		updateHistogram();

		return visualisationPanel;
	}

	private void updateHistogram() throws MicroarrayException, IOException {

		updateChipBox();

		int barCount = ((SpinnerNumberModel) barCountSpinner.getModel()).getNumber().intValue();
		String expression = ((Variable) chipBox.getSelectedItem()).getExpression();
		FloatArrayList histogram = getHistogram(barCount, expression);

		if (histogram != null) {
			CategoryDataset dataset = toCategoryDataset(histogram);

			CategoryAxis categoryAxis = new CategoryAxis("value");
			ValueAxis valueAxis = new NumberAxis("count");

			BarRenderer renderer = new BarRenderer();
			ItemLabelPosition position1 = new ItemLabelPosition(ItemLabelAnchor.OUTSIDE12, TextAnchor.BOTTOM_CENTER);
			renderer.setBasePositiveItemLabelPosition(position1);
			ItemLabelPosition position2 = new ItemLabelPosition(ItemLabelAnchor.OUTSIDE6, TextAnchor.TOP_CENTER);
			renderer.setBaseNegativeItemLabelPosition(position2);

			CategoryPlot plot = new CategoryPlot(dataset, categoryAxis, valueAxis, renderer);
			plot.setOrientation(PlotOrientation.VERTICAL);
			JFreeChart chart = new JFreeChart(data.getName(), JFreeChart.DEFAULT_TITLE_FONT, plot, false);

			visualisationPanel.removeAll();
			visualisationPanel.add(makePanel(chart), BorderLayout.CENTER);
			visualisationPanel.validate();

		} else {
			throw new IllegalArgumentException("histogram not supported for " + data.getName());
		}

	}

	@Override
	public JPanel getParameterPanel() {
		if (paramPanel == null) {
			paramPanel = new JPanel();
			paramPanel.setLayout(new GridBagLayout());
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

			chipBox = new JComboBox();

			updateButton = new JButton("Draw");
			updateButton.addActionListener(new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					try {
						// Histogram.this.updateHistogram();

						List<Variable> variables = new ArrayList<Variable>();
						variables.add((Variable) chipBox.getSelectedItem());

						application.setVisualisationMethod(new VisualisationMethodChangedEvent(this, MicroarrayModule.VisualisationMethods.HISTOGRAM, variables, getFrame().getDatas(), getFrame().getType(), getFrame()));

					} catch (Exception exp) {
						application.reportException(exp);
					}
				}
			});

			GridBagConstraints c = new GridBagConstraints();

			c.gridy = 0;
			c.gridx = 0;
			c.insets.set(5, 10, 5, 10);
			c.anchor = GridBagConstraints.NORTHWEST;
			c.weightx = 1.0;
			c.fill = GridBagConstraints.HORIZONTAL;
			paramPanel.add(new JLabel("Chip to visualise: "), c);
			c.gridy++;
			paramPanel.add(chipBox, c);
			c.gridy++;
			paramPanel.add(new JLabel("Number of bins: "), c);
			c.gridy++;
			paramPanel.add(barCountSpinner, c);
			c.gridy++;
			paramPanel.add(updateButton, c);
			c.gridy++;
			c.fill = GridBagConstraints.BOTH;
			c.weighty = 1.0;
			paramPanel.add(new JPanel(), c);
		}
		return paramPanel;
	}

	private void updateChipBox() {
		if (paramPanel == null) {
			throw new IllegalStateException("must call getParameterPanel first");
		}

		Visualisation.fillCompoBox(chipBox, this.getVariablesFor(data));

		List<Variable> variables = getFrame().getVariables();
		if (variables != null && variables.size() > 0) {
			chipBox.setSelectedItem(variables.get(0));
		}
	}

	private CategoryDataset toCategoryDataset(FloatArrayList values) {
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		for (int i = 0; i < values.size(); i++) {
			dataset.addValue(values.get(i), "data", "" + i);
		}
		return dataset;
	}

	private FloatArrayList getHistogram(int histogramSteps, String expression) throws MicroarrayException, IOException {

		// get data
		Iterable<Float> intensities = data.queryFeatures(expression).asFloats();
		LinkedList<Float> values = new LinkedList<Float>();
		for (float intensity : intensities) {
			values.add(intensity);
		}

		// sort it, so we don't have to search for value intervals
		Collections.sort(values);

		// filter out NaN's
		while (values.get(values.size() - 1).isNaN()) {
			values.remove(values.size() - 1);
		}

		if (values.size() < histogramSteps) {
			return new FloatArrayList(values); // can't make histogram, return
												// plain values
		}

		// determine step size
		float min = values.get(0);
		float max = values.get(values.size() - 1);
		float stepSize = (max - min) / ((float) histogramSteps);

		// initialise
		float[] histogram = new float[histogramSteps];
		int valueIndex = 0;
		float roof = min + stepSize;

		// step through categories, counting matching values as we go
		for (int step = 0; step < histogram.length; step++) {
			while (valueIndex < values.size() && values.get(valueIndex) <= roof) {
				histogram[step]++;
				valueIndex++;
			}
			roof += stepSize;
		}

		// add to last category what was left out
		histogram[histogram.length - 1] += (values.size() - valueIndex);

		return new FloatArrayList(histogram);
	}
}
