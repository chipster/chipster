package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.OverlayLayout;

import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.selection.RowChoiceEvent;
import fi.csc.microarray.client.selection.RowSelectionManager;
import fi.csc.microarray.client.visualisation.AnnotateListPanel;
import fi.csc.microarray.client.visualisation.ChipVisualisation;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.databeans.DataBean;

public class Scatterplot extends ChipVisualisation implements ActionListener, MouseListener, MouseMotionListener, PropertyChangeListener {

	public Scatterplot(VisualisationFrame frame) {
		super(frame);
	}

	protected ChartPanel chartPanel;

	protected JPanel paramPanel;
	protected JPanel settingsPanel;
	protected JPanel overlayPanel;
	protected AnnotateListPanel list;

	protected JComboBox xBox;
	protected JComboBox yBox;

	protected DataBean data;

	protected JButton useButton;

	@Override
	public JPanel getParameterPanel() {
		if (paramPanel == null) {
			paramPanel = new JPanel();
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);
			paramPanel.setLayout(new BorderLayout());

			JPanel settings = this.createSettingsPanel();
			list = new AnnotateListPanel();

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Settings", settings);
			tabPane.addTab("Selected", list);

			paramPanel.add(tabPane, BorderLayout.CENTER);
		}
		return paramPanel;
	}

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

		Visualisation.fillCompoBox(xBox, this.getVariablesFor(data));
		Visualisation.fillCompoBox(yBox, this.getVariablesFor(data));
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

			application.setVisualisationMethod(new VisualisationMethodChangedEvent(this, VisualisationMethod.SCATTERPLOT, vars, getFrame().getDatas(), getFrame().getType(), getFrame()));
		}
	}

	public class DataItem2D {
		private Rectangle2D bounds;
		private String name;
		private Integer index;
		private Integer series;

		public DataItem2D(Rectangle2D bounds, String name, int index, int series) {
			this.bounds = bounds;
			this.name = name;
			this.index = index;
			this.series = series;
		}

		public Rectangle2D getBounds() {
			return bounds;
		}

		public String getName() {
			return name;
		}

		public int getIndex() {
			return index;
		}

		public int getSeries() {
			return series;
		}


		public void setBounds(Rectangle rect) {
			bounds = rect;
		}

		@Override
		public boolean equals(Object other) {
			if (other instanceof DataItem2D) {
				DataItem2D otherData = (DataItem2D)other;
				return otherData.getIndex() == this.getIndex() && otherData.getSeries() == this.getSeries();
				
			} else {			
				return false;
			}
		}

		// The hashSet uses also hashCode to check equality of the objects and
		// thus overriding just the equals method isn't enought the get rid of
		// duplicates
		@Override
		public int hashCode() {
			return index.hashCode();
		}
	}

	protected List<DataItem2D> allItems = new LinkedList<DataItem2D>();
	protected Set<DataItem2D> selectedItems = new HashSet<DataItem2D>();

	protected TransparentPanel transparentPanel;

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		this.data = data;

		// Hide the panel if it exists already ( setting change) to avoid
		// repaints when the
		// allItems aren't collected yet

		if (chartPanel != null) {
			chartPanel.setVisible(false);
		}

		allItems.clear();

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

		Variable xVar = (Variable) xBox.getSelectedItem();
		Variable yVar = (Variable) yBox.getSelectedItem();

		Iterable<Float> xValues = data.queryFeatures(xVar.getExpression()).asFloats();
		Iterable<Float> yValues = data.queryFeatures(yVar.getExpression()).asFloats();

		int i = 0;
		for (String name : data.queryFeatures("/identifier").asStrings()) {
			allItems.add(new DataItem2D(null, name, i++, 0));
		}

		PlotDescription description = new PlotDescription(data.getName(), xVar.getName(), yVar.getName());

		XYSeries series = toXYSeries(xValues, yValues);
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(series);
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

	protected class TransparentPanel extends JPanel {

		private Rectangle2D area;

		public TransparentPanel() {
			this.setOpaque(false);
		}

		@Override
		public void paintComponent(Graphics g) {
			super.paintComponents(g);

			if (area != null) {
				Graphics2D g2d = (Graphics2D) g;
				g2d.setColor(Color.DARK_GRAY);
				g2d.setStroke(VisualConstants.dashLine);

				g2d.draw(area);
			}
		}

		protected void setArea(Rectangle area) {
			this.area = area;
		}
	}

	/**
	 * If x is null, enumeration from 0 is used.
	 * 
	 * @param x
	 * @param y
	 * @return
	 */
	protected XYSeries toXYSeries(Iterable<Float> xValues, Iterable<Float> yValues) {

		XYSeries series = new XYSeries("", false); // autosort=false, autosort would mess up selection
		Iterator<Float> xIterator = xValues != null ? xValues.iterator() : null;
		int i = 0;

		for (Float y : yValues) {
			if (xIterator != null) {
				series.add(xIterator.next(), y);
			} else {
				series.add(i, y);
			}
			i++;
		}
		return series;
	}

	public void mouseClicked(MouseEvent e) {
		if (!e.isControlDown()) {
			selectedItems.clear();
		}
		Rectangle2D r = null;
		for (DataItem2D item : allItems) {
			r = item.getBounds();
			if (r != null && item != null && r.contains(e.getPoint()) && !selectedItems.contains(item)) {
				selectedItems.add(item);
			}
		}

		this.list.setSelectedListContentAsDataItems(selectedItems, this, true, data);
	}

	public void mouseEntered(MouseEvent e) {
	}

	public void mouseExited(MouseEvent e) {
	}

	private Point startCoords;
	private boolean isDragged = false;

	public void mousePressed(MouseEvent e) {
		startCoords = e.getPoint();
		isDragged = false;
	}

	public void mouseReleased(MouseEvent e) {

		// Hide the selection rectangle
		transparentPanel.setArea(null);
		transparentPanel.repaint();

		if (isDragged) {
			if (!e.isControlDown()) {
				selectedItems.clear();
			}
			Rectangle2D r = null;
			int x = (int) startCoords.getX() < e.getX() ? (int) startCoords.getX() : e.getX();
			int y = (int) startCoords.getY() < e.getY() ? (int) startCoords.getY() : e.getY();
			int w = Math.abs(e.getX() - (int) startCoords.getX());
			int h = Math.abs(e.getY() - (int) startCoords.getY());

			Rectangle2D selectedArea = new Rectangle(x, y, w, h);

			for (DataItem2D item : allItems) {
				r = item.getBounds();		

				if (r != null && selectedArea.intersects(r)) {
					selectedItems.add(item);
				}
			}

			this.list.setSelectedListContentAsDataItems(selectedItems, this, true, data);
		}

		// Draw the selection frames for the dataItems
		chartPanel.getChart().fireChartChanged();
	}

	public void mouseDragged(MouseEvent e) {
		isDragged = true;

		int x = (int) startCoords.getX() < e.getX() ? (int) startCoords.getX() : e.getX();
		int y = (int) startCoords.getY() < e.getY() ? (int) startCoords.getY() : e.getY();
		int w = Math.abs(e.getX() - (int) startCoords.getX());
		int h = Math.abs(e.getY() - (int) startCoords.getY());

		transparentPanel.setArea(new Rectangle(x, y, w, h));

		transparentPanel.repaint();
	}

	public void mouseMoved(MouseEvent e) {
	}

	public void propertyChange(PropertyChangeEvent evt) {
		if (evt instanceof RowChoiceEvent && evt.getSource() != this && ((RowChoiceEvent) evt).getData() == data) {

			updateSelectionsFromApplication(false);
			chartPanel.repaint();
		}
	}

	protected void updateSelectionsFromApplication(boolean dispatchEvent) {
		RowSelectionManager manager = application.getSelectionManager().getRowSelectionManager(data);

		selectedItems.clear();
		for (Integer index : manager.getSelectedRows()) {
			DataItem2D point = allItems.get(index);
			if (point != null) {
				selectedItems.add(allItems.get(index));
			}
		}

		list.setSelectedListContentAsDataItems(selectedItems, this, dispatchEvent, data);
	}
}
