package fi.csc.microarray.client.visualisation;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.geom.Rectangle2D;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

/**
 * ChartPanel without internal scaling to be used inside JScrollPane.
 * 
 * @author klemela
 */
public class NonScalableChartPanel extends ChartPanel {
    
    public NonScalableChartPanel() {
        super(null);
    }
    
	public NonScalableChartPanel(JFreeChart chart) {
		super(chart);
	}

	@Override
	public void paintComponent(Graphics g) {

		Dimension size = getSize();
		Insets insets = getInsets();
		Rectangle2D available = new Rectangle2D.Double(insets.left, insets.top,
				size.getWidth() - insets.left - insets.right,
				size.getHeight() - insets.top - insets.bottom);
		
		this.setMinimumDrawWidth((int)available.getWidth() - 1);
		this.setMinimumDrawHeight((int)available.getHeight() - 1);
		this.setMaximumDrawWidth((int)available.getWidth() + 1);
		this.setMaximumDrawHeight((int)available.getHeight() + 1);

		super.paintComponent(g);
	}

}
