package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Shape;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import javax.swing.UIManager;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.PlotState;
import org.jfree.data.general.DatasetChangeEvent;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.EmptyTrack;

/**
 * <p>The main visual component for Genome browser. GenomePlot is the core of the view layer and runs inside Swing event
 * dispatch thread. The plot is constructed out of {@link View} components. Compatible with JFreeChart visualization library.</p>
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * @see View
 */
public class GenomePlot extends Plot {

	private List<View> views = new LinkedList<View>();
	private View dataView = null;
	private OverviewHorizontalView overviewView = null;
	private ReadScale readScale = ReadScale.AUTO;
    public ChartPanel chartPanel;
    
    private boolean showFullHeight = false;
	private Rectangle dirtyArea;
	
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

	public GenomePlot(TooltipAugmentedChartPanel panel, boolean horizontal) throws FileNotFoundException, MalformedURLException {
		
	    // set chart panel
	    this.chartPanel = panel;
	    this.chartPanel.setLayout(null);
	    this.dirtyArea = panel.getBounds();
	    
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
			public void regionChanged(Region bpRegion) {
				overviewView.highlight = bpRegion;
				overviewView.setBpRegion(new RegionDouble(0.0, 250*1000*1000.0, bpRegion.start.chr), false);
			}
		});
		
		overviewView.addOverviewRegionListener(new RegionListener() {

			@Override
			public void regionChanged(Region bpRegion) {
				dataView.setBpRegion(new RegionDouble(bpRegion), false);
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
		overviewView.setBpRegion(new RegionDouble(0d, chromosomeSizeBp, new Chromosome(chromosome)), false);
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
		RegionDouble bpCoordRegion = new RegionDouble(
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

		this.setBackgroundPaint(UIManager.getColor("Panel.background"));
		g2.setClip(this.dirtyArea);
		drawBackground(g2, this.dirtyArea); // clear everything that was drawn before 

		Shape savedClip = g2.getClip();
		g2.clip(area);

		Rectangle viewArea = (Rectangle) area.getBounds().clone();				
		this.dirtyArea = viewArea; 
		
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
				            sumStaticViewHeights()) / countStaticViews() + 1;
				}

				g2.setClip(viewArea);
				view.drawView(g2, false);
				this.dirtyArea = viewArea.union(this.dirtyArea);
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

				if (view instanceof HorizontalView) {
					viewArea.grow(-view.margin, 0);
				}

				g2.setClip(savedClip);
				g2.drawLine(viewArea.x - 1, 0, viewArea.x - 1, viewArea.height);

				g2.setClip(viewArea);
				view.drawView(g2, false);
				this.dirtyArea = viewArea.union(this.dirtyArea);
			}
		}
		g2.setClip(savedClip);
	}
	
	public int getHeightTotal() {
		int total = 0;
		for (View view : views) {
			total += view.getTrackHeightTotal();
		}
	
		return total;
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

	public void clean() {
		overviewView.clean();
		dataView.clean();
	}
}
