package fi.csc.microarray.client.visualisation;

import java.awt.Color;
import java.awt.Dimension;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JPanel;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel.SelectionChangeListener;
import fi.csc.microarray.databeans.DataBean;

public abstract class Visualisation {

	public static final Dimension PARAMETER_SIZE = new Dimension(150, 500);

	private VisualisationFrame frame;

	public abstract JComponent getVisualisation(DataBean data) throws Exception;

	// public to be able to use this also in implementations of this class and
	// in their inner classes
	public final ClientApplication application = Session.getSession().getApplication();

	public abstract boolean canVisualise(DataBean bean) throws MicroarrayException;

	public Visualisation(VisualisationFrame frame) {
		this.frame = frame;
	}

	public VisualisationFrame getFrame() {
		return frame;
	}

	public JPanel getParameterPanel() {
		return null;
	}

	/**
	 * Empty method for visualisations with several datasets. Throws
	 * NotImplementedException if this is called for the visualisation that
	 * doesn't override this method.
	 */
	public JComponent getVisualisation(List<DataBean> data) throws Exception {
		throw new NotImplementedException();
	}

	protected JComponent getDefaultVisualisation() {
		JPanel panel = new JPanel();
		panel.setBackground(Color.WHITE);
		return panel;
	}

	public static void fillCompoBox(JComboBox box, Object[] content) {
		box.removeAllItems();
		for (Object o : content) {
			box.addItem(o);
		}
	}

	/**
	 * Something in a DataBean that can be used to draw a certain kind of
	 * visualisation.
	 */
	public static class Variable {
		/**
		 * To be shown to user.
		 */
		private String name;
		/**
		 * To be be used when requesting a Feature from a DataBean.
		 */
		private String expression;

		public Variable(String name, String expression) {
			this.name = name;
			this.expression = expression;
		}

		public String getExpression() {
			return expression;
		}

		public String getName() {
			return name;
		}

		@Override
		public String toString() {
			return name;
		}

		/*
		 * Equals is overridden to recognise variables after visualisation
		 * redraw, when new variable objects area created.
		 * 
		 * @see java.lang.Object#equals(java.lang.Object)
		 */
		@Override
		public boolean equals(Object obj) {
			return obj instanceof Variable && ((Variable) obj).name.equals(this.name) && ((Variable) obj).expression.equals(this.expression);
		}
	}
	
	public Variable[] getVariablesFor(DataBean dataBean) {
		// return empty variable list
		LinkedList<Variable> vars = new LinkedList<Variable>();
		return vars.toArray(new Variable[0]);
	}

	public static ChartPanel makePanel(JFreeChart chart) {

		ChartPanel panel = new ChartPanel(chart);
		chart.getTitle().setFont(VisualConstants.VISUALISATION_TITLE_FONT);
		return panel;
	}
	
	protected JPanel makeSelectablePanel(JFreeChart chart, SelectionChangeListener selectionListener) {

		JPanel panel = new SelectableChartPanel(chart, selectionListener);
		chart.getTitle().setFont(VisualConstants.VISUALISATION_TITLE_FONT);
		return panel;
	}

	public static class PlotDescription {
		public PlotDescription(String plotTitle, String xTitle, String yTitle) {
			this.plotTitle = plotTitle;
			this.xTitle = xTitle;
			this.yTitle = yTitle;
		}

		public String plotTitle;
		public String xTitle;
		public String yTitle;
	}

	public boolean canVisualise(List<DataBean> beans) throws MicroarrayException {
		throw new NotImplementedException();
	}

	/**
	 * Visualisations that uses a dataset must override this method and return
	 * true. It will be interpreted as a commitment to override also single
	 * dataset versions of methods canVisualise and getVisualisation.
	 */

	public boolean isForSingleData() {
		return true;
	}

	/**
	 * Visualisations that use multiple datasets must override this method and
	 * return true. It will be interpreted as a commitment to override also
	 * multiple dataset versions of methods canVisualise and getVisualisation.
	 */
	public boolean isForMultipleDatas() {
		return false;
	}

	/**
	 * The references to the visualisation panel and PropetyChangeListeners will
	 * be removed automatically in VisualisationFrame.removeVisualisation().
	 * This method must be overridden if any other references are created in the
	 * visualisation. All these references have to be removed to allow garbage
	 * collection to clean visualisation.
	 */
	public void removeVisualisation() {
	}
}
