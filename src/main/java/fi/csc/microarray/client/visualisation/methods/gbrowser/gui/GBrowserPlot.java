package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Rectangle;
import java.awt.Shape;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import org.jfree.chart.plot.Plot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.chart.plot.PlotState;
import org.jfree.data.general.DatasetChangeEvent;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;

/**
 * <p>The main visual component for Genome browser. GenomePlot is the core of the view layer and runs inside Swing event
 * dispatch thread. The plot is constructed out of {@link GBrowserView} components. Compatible with JFreeChart visualization library.</p>
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * @see GBrowserViewdrawView
 */
public class GBrowserPlot extends Plot implements LayoutContainer {

	private List<GBrowserView> views = new LinkedList<GBrowserView>();
	private GBrowserView dataView = null;
	private OverviewHorizontalView overviewView = null;
	private ReadScale readScale = ReadScale.AUTO;
    public TooltipAugmentedChartPanel chartPanel;
	
	/**
	 * Scale for visualising reads as profiles, gel etc.
	 */
    public enum ReadScale {
        XS("10", 10),
        SMALL("50", 50),
        MEDIUM("100", 100),
        MEDIUM2("250", 250),
        LARGE("500", 500),        
        XL("1000", 1000),
        XXL("5000", 5000),
        XXXL("10000", 10000),
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

	public GBrowserPlot(TooltipAugmentedChartPanel panel, boolean horizontal) {
		
	    // set chart panel
	    this.chartPanel = panel;
	    this.chartPanel.setLayout(null);
	    
		// add overview view
		this.overviewView = new OverviewHorizontalView(this);
		this.overviewView.margin = 0;
		this.views.add(overviewView);

		// add horizontal or circular data view
		if (horizontal) {
			
			this.dataView = new HorizontalView(this, true, true, false);

		} else {
			
			this.dataView = new CircularView(this, true, true, false);
			this.dataView.margin = 20;
			//this.dataView.addTrack(new EmptyTrack(dataView, 30));
		}

		this.views.add(dataView);
		chartPanel.addTooltipRequestProcessor(dataView);

		dataView.addRegionListener(new RegionListener() {
			public void regionChanged(Region bpRegion) {
				overviewView.highlight = bpRegion;
				overviewView.setBpRegion(new RegionDouble(0.0, 250*1000*1000.0, bpRegion.start.chr), false);
			}
		});
		
		//Mouse click on the overview shows that region in data view
		overviewView.addOverviewRegionListener(new RegionListener() {

			@Override
			public void regionChanged(Region bpRegion) {
				dataView.setBpRegion(new RegionDouble(bpRegion), false);
			}		
		});
		
		
	}

	public GBrowserView getDataView() {
		return dataView;
	}

	public GBrowserView getOverviewView() {
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
		return "GenomeBrowser";
	}

	/**
	 * Draws the plot on a Java2D graphics device (such as the screen or a printer).
	 * 
	 * @param g2
	 *            the graphics device.
	 * @param plotArea
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
	public void draw(java.awt.Graphics2D g2, java.awt.geom.Rectangle2D plotArea, java.awt.geom.Point2D anchor, PlotState parentState, PlotRenderingInfo info) {

		if (info != null) {
			info.setPlotArea(plotArea);
			info.setDataArea(plotArea);
		}		
		
		Shape savedClip = g2.getClip();
		
		Rectangle viewArea = (Rectangle) plotArea.getBounds().clone();
		
		//Decide how to divide space between ScrollGroup, TrackGroup and Tracks
		LayoutTool.doLayout(this, (int) plotArea.getBounds().getHeight());
		
		//Draw all views
		for (int i = 0; i < views.size(); i++) {
			GBrowserView view = views.get(i);
						
			viewArea.height = view.getHeight();
			
			g2.setClip(viewArea);
			view.draw(g2, plotArea.getBounds(), viewArea);
			
			viewArea.y += viewArea.height;
		}
		
		//Update ScrollGroups' scroll bars locations and values. Height of content is known only after it is drawn.
		chartPanel.setScrollGroupOrder(getScrollGroups(), (int) plotArea.getHeight());		
		
		g2.setClip(savedClip);
	}

	private Collection<ScrollGroup> getScrollGroups() {
		
		List<ScrollGroup> groups = new LinkedList<ScrollGroup>();
		for ( GBrowserView view : views) {
			groups.addAll(view.getScrollGroups());
		}
		return groups;
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
		if (!(obj instanceof GBrowserPlot)) {
			return false;
		}

		// can't find any difference...
		return true;
	}

	public void redraw() {
		this.datasetChanged(new DatasetChangeEvent(this, null));		
	}

	public Collection<GBrowserView> getViews() {
		return views;
	}
	
    public ReadScale getReadScale() {
        return readScale;
    }

    public void setReadScale(ReadScale readScale) {
        this.readScale = readScale;
        this.dataView.redraw();
    }
    
    public void setFullLayoutMode(boolean enabled) {
    	
    	LayoutTool.setFullLayoutMode(this, enabled);
    
    	redraw();
    }

	public void clean() {		
		chartPanel.clean();
		overviewView.clean();
		dataView.clean();
		dataView = null;
	}

	@Override
	public Collection<? extends LayoutComponent> getLayoutComponents() {
		return views;
	}

	public void initializeTracks() {
		for (GBrowserView view : views) {
			for (Track track : view.getTracks()) {
				track.initializeListener();
			}			
		}		
	}

	public void updateData() {
		for (GBrowserView view : views) {
			view.fireAreaRequests();			
		}
	}
}
