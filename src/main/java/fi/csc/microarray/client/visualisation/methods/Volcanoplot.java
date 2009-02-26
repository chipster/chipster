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

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel.SelectionChangeListener;
import fi.csc.microarray.databeans.DataBean;

public class Volcanoplot extends Scatterplot implements ActionListener, PropertyChangeListener, SelectionChangeListener {

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
		
		refreshAxisBoxes(data);
		
		List<Variable> vars = getFrame().getVariables();
		
		//If this a redraw from the settings panel, use asked columns
		if (vars != null && vars.size() == 2) {			
			xBox.setSelectedItem(vars.get(0));
			yBox.setSelectedItem(vars.get(1));
		}

		xVar = (Variable) xBox.getSelectedItem();
		yVar = (Variable) yBox.getSelectedItem();

		PlotDescription description = new PlotDescription(data.getName(), "fold change", "-log(p)");

		NumberAxis domainAxis = new NumberAxis(description.xTitle);
		NumberAxis rangeAxis = new NumberAxis(description.yTitle);
		rangeAxis.setRange(new Range(0, 16));
		
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

		JFreeChart chart = new JFreeChart(description.plotTitle, plot);

		application.addPropertyChangeListener(this);

		selectableChartPanel = new SelectableChartPanel(chart, this); 
		return selectableChartPanel;
	}
	
	private Iterator<Float> getXValueIterator() throws MicroarrayException{

		return data.queryFeatures(xVar.getExpression()).asFloats().iterator();		
	}

	private Iterator<Float> getYValueIterator() throws MicroarrayException{

		// "/column/" part of the query comes from the getExpression function		
		return  data.queryFeatures(Y_AXIS_CONVERTER.replace(
				COLUMN_MASK, yVar.getExpression())).asFloats().iterator();
	}
	
	protected void updateXYSerieses() throws MicroarrayException {
		

		Iterator<Float> xValues = getXValueIterator();
		Iterator<Float> yValues = getYValueIterator();
		
		XYSeries greenSeries = new XYSeries(""); 
		XYSeries blackSeries = new XYSeries(""); 
		XYSeries redSeries = new XYSeries(""); 
		XYSeries selectedSeries = new XYSeries("");
		
		int row = 0;
		while( xValues.hasNext() && yValues.hasNext()) {
			
			float x = xValues.next();
			float y = yValues.next();
			boolean overYThreshold = y >= -Math.log(0.05);
			boolean overXThreshold = Math.abs(x) >= 1f;
			
			if( selectedIds.contains(row)){
				selectedSeries.add(new XYDataItem(x, y));
			} else {

				if (overYThreshold && overXThreshold) {
					if(x < 0){
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
	
public void selectionChanged(Rectangle.Double newSelection) {
		
		if(newSelection == null){
			selectedIds.clear();
		} else {
		
			Iterator<Float> xValues;
			Iterator<Float> yValues;
			
			try {								
				
				xValues = getXValueIterator();
				yValues = getYValueIterator();

				for (int i = 0;	xValues.hasNext() && yValues.hasNext();	i++){			

					if(newSelection.contains(new Point.Double(xValues.next(), yValues.next()))){

						//Contains method should work with Integers as it uses equals to compare objects. 
						//Usage of hash can be still a problem, as VM pools integer objects only for 
						//integers between -256 and 256 or something like that.
						if(selectedIds.contains(i)){
							//Remove from selection if selected twice
							selectedIds.remove(i);
						} else {
							selectedIds.add(i);
						}
					}
				}		
			} catch (MicroarrayException e) {
				application.reportException(e);
			}
		}
		
		this.list.setSelectedRows(selectedIds, this, true, data);
		
		try {
			updateXYSerieses();
		} catch (MicroarrayException e) {
			application.reportException(e);
		}
	}
}
