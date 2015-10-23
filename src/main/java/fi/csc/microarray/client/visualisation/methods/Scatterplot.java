package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.Ellipse2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import fi.csc.microarray.client.selection.SelectionEvent;
import fi.csc.microarray.client.selection.IntegratedSelectionManager;
import fi.csc.microarray.client.visualisation.SelectionList;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel.SelectionChangeListener;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class Scatterplot extends ChipVisualisation 
implements ActionListener, PropertyChangeListener, SelectionChangeListener {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	protected SelectableChartPanel selectableChartPanel;

	protected JPanel paramPanel;
	protected SelectionList list;

	protected JComboBox<Variable> xBox;
	protected JComboBox<Variable> yBox;
	
	protected Variable xVar;
	protected Variable yVar;
	
	protected XYPlot plot;

	protected DataBean data;

	protected JButton useButton;

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
		settingsPanel.add(new JLabel("X-axis: "), c);
		c.gridy++;
		settingsPanel.add(xBox, c);
		c.gridy++;
		settingsPanel.add(new JLabel("Y-axis: "), c);
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

		Visualisation.fillComboBox(xBox, this.getVariablesFor(data));
		Visualisation.fillComboBox(yBox, this.getVariablesFor(data));
	}

	/**
	 * A method defined by the ActionListener interface. Allows this panel to
	 * listen to actions on its components.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == useButton) {
			List<Variable> vars = new ArrayList<Variable>();
			vars.add((Variable) xBox.getSelectedItem());
			vars.add((Variable) yBox.getSelectedItem());

			application.setVisualisationMethod(new VisualisationMethodChangedEvent(this, MicroarrayModule.VisualisationMethods.SCATTERPLOT, vars, getFrame().getDatas(), getFrame().getType(), getFrame()));
		}
	}
	
	protected Set<Integer> selectedIndexes = new HashSet<Integer>();

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		this.data = data;

		refreshAxisBoxes(data);

		List<Variable> vars = getFrame().getVariables();
		if (vars == null || vars.size() < 2) {
			if (xBox.getItemCount() >= 1) {
				xBox.setSelectedIndex(0);
			}
			if (yBox.getItemCount() >= 2) {
				yBox.setSelectedIndex(1);
			} else {
				yBox.setSelectedIndex(0);
			}
		} else {
			xBox.setSelectedItem(vars.get(0));
			yBox.setSelectedItem(vars.get(1));
		}

		xVar = (Variable) xBox.getSelectedItem();
		yVar = (Variable) yBox.getSelectedItem();

		PlotDescription description = new PlotDescription(data.getName(), xVar.getName(), yVar.getName());

		NumberAxis domainAxis = new NumberAxis(description.xTitle);
		domainAxis.setAutoRangeIncludesZero(false);
		NumberAxis rangeAxis = new NumberAxis(description.yTitle);
		rangeAxis.setAutoRangeIncludesZero(false);
		
		XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();		
		renderer.setLinesVisible(false);
		renderer.setShapesVisible(true);
		renderer.setShape(new Ellipse2D.Float(-2, -2, 4, 4));
		renderer.setSeriesPaint(1, Color.black);
		
		plot = new XYPlot(new XYSeriesCollection(), domainAxis, rangeAxis, renderer);
		
		this.updateSelectionsFromApplication(false);		//Calls also updateXYSerieses();
		
		JFreeChart chart = new JFreeChart(description.plotTitle, plot);

		application.addClientEventListener(this);

		selectableChartPanel = new SelectableChartPanel(chart, this); 
		return selectableChartPanel;
	}


	protected void updateXYSerieses() throws MicroarrayException {

		Iterable<Float> xValues = data.queryFeatures(xVar.getExpression()).asFloats();
		Iterable<Float> yValues = data.queryFeatures(yVar.getExpression()).asFloats();
		
		XYSeries series = new XYSeries(""); 
		XYSeries selectionSeries = new XYSeries("");
		
		Iterator<Float> xIterator = xValues != null ? xValues.iterator() : null;
		int i = 0;

		for (Float y : yValues) {
			if(selectedIndexes.contains(i)){
				if (xIterator != null) {
					selectionSeries.add(xIterator.next(), y);
				} else {
					selectionSeries.add(i, y);
				}
			} else {

				if (xIterator != null) {
					series.add(xIterator.next(), y);
				} else {
					series.add(i, y);
				}
			}
			i++;
		}
		
			
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
		dataset.addSeries(selectionSeries);

		plot.setDataset(dataset);

	}

	public void propertyChange(PropertyChangeEvent evt) {
		if (evt instanceof SelectionEvent && evt.getSource() != this && ((SelectionEvent) evt).getData() == data) {

			updateSelectionsFromApplication(false);
		}
	}

	protected void updateSelectionsFromApplication(boolean dispatchEvent) {
		IntegratedSelectionManager manager = application.getSelectionManager().getSelectionManager(data);

		selectedIndexes.clear();
		for (int i : manager.getSelectionAsRows()){
			selectedIndexes.add(i);
		}

		list.setSelectedRows(selectedIndexes, this, dispatchEvent, data);
		
		try {
			updateXYSerieses();
		} catch (MicroarrayException e) {
			application.reportException(e);
		}
	}

	public void selectionChanged(Rectangle.Double newSelection) {
		
		if(newSelection == null){
			selectedIndexes.clear();
		} else {
		
			Iterator<Float> xValues;
			Iterator<Float> yValues;
			
			try {								
				
				xValues = data.queryFeatures(xVar.getExpression()).asFloats().iterator();
				yValues = data.queryFeatures(yVar.getExpression()).asFloats().iterator();

				for (int i = 0;	xValues.hasNext() && yValues.hasNext();	i++){			

					double x = xValues.next();
					double y = yValues.next();				
					
					if(newSelection.contains(new Point.Double(x, y))){

						//Contains method should work with Intgers as it uses equals to compare objects. 
						//Usage of hash inside equals can be still a problem, as VM pools integer 
						//objects only for integers between -256 and 256 or something like that.
						if(selectedIndexes.contains(i)){
							//Remove from selection if selected twice
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
