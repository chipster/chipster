package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Shape;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.net.MalformedURLException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import org.jfree.chart.ChartMouseEvent;
import org.jfree.chart.ChartMouseListener;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.PlotState;
import org.jfree.data.general.DatasetChangeEvent;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.EmptyTrack;

/**
 * The main visual component for Genome Browser. Compatible with JFreeChart. 
 * 
 * @author Petri Klemelä, Aleksi Kallio
 */
public class GenomePlot extends Plot implements ChartMouseListener, Cloneable, Serializable {

	private List<View> views = new LinkedList<View>();
	private View dataView = null;
	private OverviewHorizontalView overviewView = null;
	private ReadScale readScale = ReadScale.AUTO;
    public ChartPanel chartPanel;
    
    private boolean showFullHeight = false;
	
	/**
	 * Scale for visualising reads as profiles, gel etc.
	 */
    public enum ReadScale {
        XS("0..10", 10),
        SMALL("0..50", 50),
        MEDIUM("0..100", 100),
        LARGE("0..500", 500),        
        XL("0..1000", 1000),
        XXL("0..5000", 5000),
        XXXL("0..10000", 10000),
        AUTO("Automatic", 0);
        
        private String name;
        public Integer numReads;
        
        private ReadScale(String name, Integer numReads) {
            this.name = name;
            this.numReads = numReads;
        }
        
        public String toString() {
            return name;
        }
    }    

	public GenomePlot(TooltipEnabledChartPanel panel, boolean horizontal) throws FileNotFoundException, MalformedURLException {
	    
	    // set chart panel
	    chartPanel = panel;
	    chartPanel.setLayout(null);

		// add overview view
		this.overviewView = new OverviewHorizontalView(this);
		this.overviewView.margin = 0;
		this.overviewView.setStaticHeight(true);
		this.overviewView.setHeight(25);
		this.views.add(overviewView);

		// add horizontal or circular data view
		if (horizontal) {
			
			this.dataView = new HorizontalView(this, true, true, false);

		} else {
			
			this.dataView = new CircularView(this, true, true, false);
			this.dataView.margin = 20;
			this.dataView.addTrack(new EmptyTrack(dataView, 30));
		}

		this.views.add(dataView);
		panel.addTooltipRequestProcessor(dataView);

		dataView.addRegionListener(new RegionListener() {
			public void regionChanged(BpCoordRegion bpRegion) {
				overviewView.highlight = bpRegion;
				overviewView.setBpRegion(new BpCoordRegionDouble(0.0, 250*1000*1000.0, bpRegion.start.chr), false);
			}
		});
		
		overviewView.addOverviewRegionListener(new RegionListener() {

			@Override
			public void regionChanged(BpCoordRegion bpRegion) {
				dataView.setBpRegion(new BpCoordRegionDouble(bpRegion), false);
			}		
		});
	}

	public View getDataView() {
		return dataView;
	}

	public View getOverviewView() {
		return overviewView;
	}

	/**
	 * Start genome plot by moving both overview and data views to some location.
	 * 
	 * @param chromosome chromosome identifier
	 * @param chromosomeSizeBp chromosome size in base pairs
	 * @param position which bp should be in the middle
	 * @param length number of visible base pairs (zoom)
	 */
	public void start(String chromosome, Double chromosomeSizeBp, Long position, Long length) {
		overviewView.setBpRegion(new BpCoordRegionDouble(0d, chromosomeSizeBp, new Chromosome(chromosome)), false);
		moveDataBpRegion(new Chromosome(chromosome), position, length);
	}

    /**
     * Move data view to a certain location in the genome.
     * 
     * @param moveToChr
     * @param moveToBp
     * @param length
     */
	public void moveDataBpRegion(Chromosome moveToChr, Long moveToBp, Long length) {
		BpCoordRegionDouble bpCoordRegion = new BpCoordRegionDouble(
				new Double(moveToBp - (length/2)),
				new Double(moveToBp + (length/2)), 
				moveToChr
		);
		dataView.setBpRegion(bpCoordRegion, false);
	}
	
	public void addDataRegionListener(RegionListener regionListener) {
		dataView.addRegionListener(regionListener);		
	}
	
	public String getPlotType() {
		return "GeneBrowser";
	}

	/**
	 * Draws the plot on a Java2D graphics device (such as the screen or a printer).
	 * 
	 * @param g2
	 *            the graphics device.
	 * @param area
	 *            the area within which the plot should be drawn.
	 * @param anchor
	 *            the anchor point (<code>null</code> permitted).
	 * @param parentState
	 *            the state from the parent plot, if there is one (<code>null</code> permitted.)
	 * @param info
	 *            collects info about the drawing (<code>null</code> permitted).
	 * @throws NullPointerException
	 *             if g2 or area is null.
	 */
	public void draw(java.awt.Graphics2D g2, java.awt.geom.Rectangle2D area, java.awt.geom.Point2D anchor, PlotState parentState, PlotRenderingInfo info) {

		if (info != null) {
			info.setPlotArea(area);
			info.setDataArea(area);
		}

		this.setBackgroundPaint(Color.black);

		drawBackground(g2, area);
		drawOutline(g2, area);

		Shape savedClip = g2.getClip();
		g2.clip(area);

		Rectangle viewArea = (Rectangle) area.getBounds().clone();				

		// Horizontal or vertical split
		if (true) {

			for (int i = 0; i < views.size(); i++) {
			    View view = views.get(i);
			    
				if (i > 0) {
					viewArea.y += viewArea.height;
				}
				
				if (view.hasStaticHeight()) {
				    viewArea.height = (int) (view.getHeight());
				} else {
				    viewArea.height = (int) (area.getBounds().getHeight() -
				            sumStaticViewHeights()) / countStaticViews();
				}

				g2.setClip(viewArea);
				view.drawView(g2, false);
			}

		} else {
			float[] viewWidths = new float[] { 0.05f, 0.95f };
			Rectangle lastArea = null;

			for (int i = 0; i < views.size(); i++) {
				if (lastArea != null) {
					viewArea.x = lastArea.x + lastArea.width;
					viewArea.y = lastArea.y;
					viewArea.height = lastArea.height;
				}
				g2.setColor(Color.black);
				viewArea.width = (int) (area.getBounds().getWidth() * viewWidths[i]);
				lastArea = (Rectangle) (viewArea.clone());

				View view = views.get(i);

				if (view instanceof VerticalView) {
					viewArea.grow(0, -view.margin);
				} else if (view instanceof HorizontalView) {
					viewArea.grow(-view.margin, 0);
				}

				g2.setClip(savedClip);
				g2.drawLine(viewArea.x - 1, 0, viewArea.x - 1, viewArea.height);

				g2.setClip(viewArea);
				view.drawView(g2, false);
			}
		}
		g2.setClip(savedClip);

		drawOutline(g2, area);
	}
	
	public int getHeightTotal() {
		int total = 0;
		for (View view : views) {
			total += view.getTrackHeightTotal();
		}
	
		return total;
	}

	/**
	 * Implements the ChartMouseListener interface. This method does nothing.
	 * 
	 * @param event
	 *            the mouse event.
	 */
	public void chartMouseMoved(ChartMouseEvent event) {
	}

	public void chartMouseClicked(ChartMouseEvent e) {
	}

	/**
	 * Tests this plot for equality with an arbitrary object. Note that the plot's dataset is NOT included in the test for equality.
	 * 
	 * @param obj
	 *            the object to test against (<code>null</code> permitted).
	 * 
	 * @return <code>true</code> or <code>false</code>.
	 */
	public boolean equals(Object obj) {
		if (obj == this) {
			return true;
		}
		if (!(obj instanceof GenomePlot)) {
			return false;
		}

		// can't find any difference...
		return true;
	}

	/**
	 * Provides serialization support.
	 * 
	 * @param stream
	 *            the output stream.
	 * 
	 * @throws IOException
	 *             if there is an I/O error.
	 * @throws NullPointerException
	 *             if stream is null.
	 */
	private void writeObject(ObjectOutputStream stream) throws IOException {
		stream.defaultWriteObject();
	}

	/**
	 * Provides serialization support.
	 * 
	 * @param stream
	 *            the input stream.
	 * 
	 * @throws IOException
	 *             if there is an I/O error.
	 * @throws ClassNotFoundException
	 *             if there is a classpath problem.
	 * @throws NullPointerException
	 *             if stream is null.
	 */
	private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
		stream.defaultReadObject();
	}

	public void redraw() {
		this.datasetChanged(new DatasetChangeEvent(this, null));
	}

	public Collection<View> getViews() {
		return views;
	}
	
    public ReadScale getReadScale() {
        return readScale;
    }

    public void setReadScale(ReadScale readScale) {
        this.readScale = readScale;
    }
    
    /**
     * Sum heights of all views in this plot that have constant heights.
     */
    private int sumStaticViewHeights() {
        int heightSum = 0;
        for (View view : views) {
            if (view.hasStaticHeight()) {
                heightSum += view.getHeight();
            }
        }
        return heightSum;
    }
    
    /**
     * Return a number of views that have static heights.
     */
    private int countStaticViews() {
        int count = 0;
        for (View view : views) {
            if (view.hasStaticHeight()) {
                count += 1;
            }
        }
        return count;
    }
    
    public boolean isFullHeight() {
    	return showFullHeight;
    }
    
    public void setFullHeight(boolean b) {
    	showFullHeight = b;
    }
}
