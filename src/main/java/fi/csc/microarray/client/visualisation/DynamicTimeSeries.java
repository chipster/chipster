package fi.csc.microarray.client.visualisation;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.DateAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.data.time.Millisecond;
import org.jfree.data.time.TimeSeries;
import org.jfree.data.time.TimeSeriesCollection;

import fi.csc.microarray.util.FloatArrayList;

public class DynamicTimeSeries {
	private JPanel panel;
	private List<TimeSeries> serieses = new LinkedList<TimeSeries>();

	public DynamicTimeSeries(List<String> seriesNames, int historyCount, String title) {
		TimeSeriesCollection dataset = new TimeSeriesCollection(); 
		for (String name : seriesNames) {
			TimeSeries ts = new TimeSeries(name, Millisecond.class);
			ts.setMaximumItemCount(historyCount);
			dataset.addSeries(ts);
			serieses.add(ts);
		}

		DateAxis domain = new DateAxis("Time");
		NumberAxis range = new NumberAxis(title);
		domain.setTickLabelsVisible(true);
		domain.setAutoRange(true);
		range.setAutoRange(true);
		XYItemRenderer renderer = new XYLineAndShapeRenderer(true, false);
		renderer.setSeriesPaint(0, Color.RED);
		renderer.setSeriesPaint(1, Color.RED);
		renderer.setStroke(new BasicStroke(3f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_BEVEL));
		XYPlot plot = new XYPlot(dataset, domain, range, renderer);
		JFreeChart chart = new JFreeChart(title, new Font("SansSerif", Font.BOLD, 23), plot, true);
		this.panel = new ChartPanel(chart);
	}

	public JPanel getPanel() {
		return panel;
	}
	
	public void addTimePoint(FloatArrayList values) {
		assert(values.size() == serieses.size());
		int i = 0;
		for (TimeSeries ts : serieses) {
			ts.add(new Millisecond(), values.get(i));
			i++;
		}		
	}
}
