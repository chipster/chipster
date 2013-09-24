package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Component;
import java.awt.Cursor;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.io.IOException;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JFileChooser;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.util.ImageExportUtils;

/**
 * <p>The main visual component for Genome browser. GenomePlot is the core of the view layer and runs inside Swing event
 * dispatch thread. The plot is constructed out of {@link GBrowserView} components. Implements also support for printing the component
 * or saving it as an image.</p>
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * @see GBrowserViewdrawView
 */
public class GBrowserPlot implements ActionListener, Printable {
	
	JPanel component = new JPanel() {
		
	};
	
	private List<GBrowserView> views = new LinkedList<GBrowserView>();
	private GBrowserView dataView = null;
	private OverviewHorizontalView overviewView = null;
	private ReadScale readScale = ReadScale.AUTO;
	
	private JMenuItem saveMenuItem;
	private JMenuItem printMenuItem;
	private GBrowser browser;
	private JFileChooser saveFileChooser;
	
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

	public GBrowserPlot(GBrowser browser, boolean horizontal, Region defaultLocation) {		
				
		this.browser = browser;
	    
		// add overview view
		this.overviewView = new OverviewHorizontalView(this);
		this.dataView = new HorizontalView(this, true, true, false);

		this.views.add(overviewView);
		this.views.add(dataView);
		
		component.setLayout(new MigLayout("flowy, fillx, gap 0! 0!, insets 0"));
		component.add(overviewView.getComponent(), "growx");					
		component.add(dataView.getComponent(), "grow");		
				
//		chartPanel.addTooltipRequestProcessor(dataView);

		dataView.addRegionListener(new RegionListener() {
			public void regionChanged(Region bpRegion) {
						
				//Change chromosome
				//This region is bigger than real chromosomes and it will be limited to actual size of the chromosome
				overviewView.setBpRegion(new RegionDouble(-50*1000*1000.0d, 300*1000*1000.0d, bpRegion.start.chr));
				overviewView.highlight = bpRegion;
			}
		});
		
		//Mouse click on the overview shows that region in data view
		overviewView.addOverviewRegionListener(new RegionListener() {

			@Override
			public void regionChanged(Region bpRegion) {
				dataView.setBpRegion(new RegionDouble(bpRegion));
			}		
		});
		
		component.setCursor(new Cursor(Cursor.HAND_CURSOR));
		
		JPopupMenu popup = new JPopupMenu();
		saveMenuItem = new JMenuItem("Save as...");
		printMenuItem = new JMenuItem("Print...");
		saveMenuItem.addActionListener(this);
		printMenuItem.addActionListener(this);
		popup.add(saveMenuItem);
		popup.add(printMenuItem);
		component.setComponentPopupMenu(popup);
		
		//Set default location to plot to avoid trouble in track initialization. 
		dataView.setBpRegion(new RegionDouble(defaultLocation));
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
		overviewView.setBpRegion(new RegionDouble(0d, chromosomeSizeBp, new Chromosome(chromosome)));
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
		dataView.setBpRegion(bpCoordRegion);
	}
	
	public void addDataRegionListener(RegionListener regionListener) {
		dataView.addRegionListener(regionListener);		
	}
	
	public String getPlotType() {
		return "GenomeBrowser";
	}

	public void redraw() {
		for (GBrowserView view : views) {
			view.updateLayout();
		}
		component.repaint();
	}

	public Collection<GBrowserView> getViews() {
		return views;
	}
	
    public ReadScale getReadScale() {
        return readScale;
    }

    public void setReadScale(ReadScale readScale) {
        this.readScale = readScale;
        this.dataView.reloadData();
    }

	public void clean() {		
		overviewView.clean();
		dataView.clean();
		dataView = null;
	}

	public void initializeDataResultListeners() {
		for (GBrowserView view : views) {
			for (TrackGroup group : view.getTrackGroups()) {
				group.initializeListener();
			}			
		}		
	}

	public void updateData() {
		for (GBrowserView view : views) {
			view.fireDataRequests();			
		}
	}

	@Override
	public void actionPerformed(ActionEvent event) {
		if (event.getSource() == saveMenuItem) {
			try {
				saveFileChooser = ImageExportUtils.saveComponent(component, saveFileChooser);
				//Remove the  visual artifacts caused by the rendering of external buffer
				this.redraw();
			} catch (IOException e) {
				browser.reportException(e);
			}
		}
		
		if (event.getSource() == printMenuItem) {

			try {
				ImageExportUtils.printComponent(this);
				//Remove the  visual artifacts caused by the rendering of external buffer
				this.redraw();
			} catch (PrinterException e) {	
				browser.reportException(e);
			}
		}
	}

	@Override
	public int print(Graphics g, PageFormat pf, int page)
			throws PrinterException {
		return ImageExportUtils.printComponent(g, pf, page, component);		
	}

	public Component getComponent() {
		return component;
	}

	public SelectionManager getSelectionManager() {
		return browser.getSelectionManager();
	}

	public GBrowser getBrowser() {
		return browser;
	}
}
