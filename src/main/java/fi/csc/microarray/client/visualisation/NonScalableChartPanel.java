package fi.csc.microarray.client.visualisation;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.geom.Rectangle2D;

import javax.swing.JViewport;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomePlot;

/**
 * ChartPanel without internal scaling to be used inside JScrollPane.
 * 
 * @author klemela
 */
public class NonScalableChartPanel extends ChartPanel {
	
	private GenomePlot genomePlot;
    
    public NonScalableChartPanel() {
        super(null);
    }
    
	public NonScalableChartPanel(JFreeChart chart) {
		super(chart);
	}

	@Override
	public void paintComponent(Graphics g) {

		//in genomeBrowser getSize() doesn't work below width of 680
		//but parent seems to know better
		
		int height =  getParent().getSize().height;
		
		/* If the genomeBrowser height is bigger than the scrollPane height, we set that bigger value
		 * to both chartpanel size and drawHeight to make the genomeBrowser draw vertically everything
		 * and still no scale it.
		 */
		if (genomePlot != null && genomePlot.isFullHeight()) {
						
			height = genomePlot.getHeight();
		}
		
		if (getParent() instanceof JViewport) {
			JViewport viewport = (JViewport)getParent();
			genomePlot.setFullHeightClip(viewport.getViewRect());
		}
				
		Dimension size = new Dimension(getParent().getSize().width, height);
		
		//Required for JScrollPane to understand that content requires scrolling
		this.setPreferredSize(size);
		
		//Required for JFreeChart to avoid scaling after the view is drawn
		this.setSize(size);
				
		//Required for JFreeChart to do the view drawing in correct resolution
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

	public void setGenomePlot(GenomePlot plot) {
		genomePlot = plot;		
	}

}
