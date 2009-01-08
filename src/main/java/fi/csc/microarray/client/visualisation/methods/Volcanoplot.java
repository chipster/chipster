package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Ellipse2D;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
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
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.Range;
import org.jfree.data.xy.XYDataItem;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.OverlayLayout;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.databeans.DataBean;

public class Volcanoplot extends Scatterplot implements ActionListener, MouseListener, MouseMotionListener, PropertyChangeListener {

	/**
	 * Analytical expression is -log(p), where p is the first available p.* column.
	 */
	private static final String COLUMN_MASK = "***";
	private static final String Y_AXIS_CONVERTER= "neg(log(***))";
	private static final String Y_AXIS_COLUMN_HEADER = "p.";
	private static final String X_AXIS_COLUMN_HEADER = "FC";

	public Volcanoplot(VisualisationFrame frame) {
		super(frame);
	}
	
	@Override
	public JPanel createSettingsPanel() {

		settingsPanel = new JPanel();
		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

		xBox = new JComboBox();
		yBox = new JComboBox();

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

		Visualisation.fillCompoBox(xBox, VisualisationUtilities.getVariablesFiltered(
				data, X_AXIS_COLUMN_HEADER, false));
		Visualisation.fillCompoBox(yBox, VisualisationUtilities.getVariablesFiltered(
				data, Y_AXIS_COLUMN_HEADER, false));
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
		
		refreshAxisBoxes(data);
		
		List<Variable> vars = getFrame().getVariables();
		
		//If this a redraw from the settings panel, use asked columns
		if (vars != null && vars.size() == 2) {			
			xBox.setSelectedItem(vars.get(0));
			yBox.setSelectedItem(vars.get(1));
		}

		Variable xVar = (Variable) xBox.getSelectedItem();
		Variable yVar = (Variable) yBox.getSelectedItem();

		// "/column/" part of the query comes from the getExpression function
		Iterator<Float> xValues = data.queryFeatures(
				xVar.getExpression()).asFloats().iterator();
		Iterator<Float> yValues = data.queryFeatures(
				Y_AXIS_CONVERTER.replace(COLUMN_MASK, yVar.getExpression())).asFloats().iterator();

		XYSeries greenSeries = new XYSeries("", false); // autosort=false, autosort would mess up selection
		XYSeries blackSeries = new XYSeries("", false); // autosort=false, autosort would mess up selection
		XYSeries redSeries = new XYSeries("", false); // autosort=false, autosort would mess up selection
		int row = 0;
		for (String name : data.queryFeatures("/identifier").asStrings()) {
			float x = xValues.next();
			float y = yValues.next();
			boolean overYThreshold = y >= -Math.log(0.05);
			boolean overXThreshold = Math.abs(x) >= 1f;
			int series;
			int index;
			
			if (overYThreshold && overXThreshold) {
				if(x < 0){
					series = 0;
					index = greenSeries.getItemCount();
					greenSeries.add(new XYDataItem(x, y));
					
				} else {			
					series = 1;
					index = redSeries.getItemCount();
					redSeries.add(new XYDataItem(x, y));
					
				}
			} else {
				series = 2;
				index = blackSeries.getItemCount();
				blackSeries.add(new XYDataItem(x, y));
			
			}
			allItems.add(new DataItem2D(null, name, row, index, series));
			row++;
		}

		PlotDescription description = new PlotDescription(data.getName(), "fold change", "-log(p)");

		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(greenSeries);
		dataset.addSeries(redSeries);
		dataset.addSeries(blackSeries);

		NumberAxis domainAxis = new NumberAxis(description.xTitle);
		NumberAxis rangeAxis = new NumberAxis(description.yTitle);
		rangeAxis.setRange(new Range(0, 16));
		
		
		XYPlot plot = new XYPlot(dataset, domainAxis, rangeAxis, null);
		JFreeChart chart = new JFreeChart(description.plotTitle, plot);
		chartPanel = makePanel(chart);
		chartPanel.addMouseListener(this);
		chartPanel.addMouseMotionListener(this);
		chartPanel.setMouseZoomable(false);
		
		// rendered depends on chartPanel for coordinate translations (a bit awkward...), so it cannot be created earlier
		XYItemRenderer renderer = new PositionRecordingRenderer(StandardXYItemRenderer.SHAPES, allItems, selectedItems, chartPanel);
		
		renderer.setSeriesPaint(0, Color.green);
		renderer.setSeriesPaint(1, Color.red);
		renderer.setSeriesPaint(2, Color.black);
		renderer.setShape(new Ellipse2D.Float(-2, -2, 4, 4));
		
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
		return isTabular && hasRows(bean) && 
			bean.queryFeatures(Y_AXIS_CONVERTER.replace(
					COLUMN_MASK, "/column/" + Y_AXIS_COLUMN_HEADER + "*")).exists() && 
			bean.queryFeatures("/column/" + X_AXIS_COLUMN_HEADER + "*").exists();
	}
	
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == useButton) {
			List<Variable> vars = new ArrayList<Variable>();
			vars.add((Variable) xBox.getSelectedItem());
			vars.add((Variable) yBox.getSelectedItem());

			application.setVisualisationMethod(new VisualisationMethodChangedEvent(
					this, VisualisationMethod.VOLCANOPLOT, vars, getFrame().getDatas(), getFrame().getType(), getFrame()));
		}
	}
}
