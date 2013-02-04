package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.CytobandDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;

/**
 * Utility class for creating predefined {@link TrackGroup} objects.  
 *
 *@author Petri Klemelä, Aleksi Kallio
 *
 */
public class TrackFactory {
	
    public static TrackGroup getGeneTrackGroup(GBrowserPlot genomePlot, DataSource annotationDataSource, TabixDataSource repeatDataSource) throws FileNotFoundException {
        
		GBrowserView dataView = genomePlot.getDataView();
		
		TrackGroup geneGroup = new GeneTrackGroup(dataView, annotationDataSource, repeatDataSource);
			    
	    return geneGroup;
	}
    
	public static TrackGroup getThinSeparatorTrackGroup(GBrowserPlot genomePlot) {
		GBrowserView view = genomePlot.getDataView();
		TrackGroup group = new TrackGroup(view);
		group.addTrack(new SeparatorTrack(view, Color.LIGHT_GRAY, 3, 0, Long.MAX_VALUE));
		return group;
	}

	public static TrackGroup getThickSeparatorTrackGroup(GBrowserPlot genomePlot) {
		GBrowserView view = genomePlot.getDataView();
		TrackGroup group = new TrackGroup(view);
		group.addTrack(new SeparatorTrack3D(view, 0, Long.MAX_VALUE, false));
		group.addTrack(new EmptyTrack(view, 2));
		group.addTrack(new SeparatorTrack3D(view, 0, Long.MAX_VALUE, true));
		
		return group;
	}
	
	public static TrackGroup getReadTrackGroup(GBrowserPlot genomePlot, DataSource userData, DataSource seqFile, String title)
	        throws FileNotFoundException, MalformedURLException {
	
		GBrowserView dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadTrackGroup readGroup = new ReadTrackGroup(dataView, userData, seqFile, title);
		readGroup.initialise();
        
        return readGroup;
	}

	public static TrackGroup getReadSummaryTrackGroup(GBrowserPlot genomePlot, DataSource userData,
	        DataSource seqFile, String title, TabixDataSource summaryDataSource)
	        throws FileNotFoundException, MalformedURLException {
	
		GBrowserView dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadSummaryTrackGroup readGroup = new ReadSummaryTrackGroup(dataView, userData, seqFile, title, summaryDataSource);
		readGroup.initialise();
        
        return readGroup;
	}

	public static TrackGroup getWigTrackGroup(GBrowserPlot plot, DataSource peakFile) {
		WIGTrack annotation = new WIGTrack(plot.getDataView(), peakFile,
		        Color.BLUE, 0, Long.MAX_VALUE);
		return new TrackGroup(annotation);
	}
	
	public static TrackGroup getPeakTrackGroup(GBrowserPlot plot, DataSource peaks) {
		GBrowserView dataView = plot.getDataView();

		PeakTrack annotation = new PeakTrack(dataView, peaks, GBrowserConstants.BED_COLOR, 0, Long.MAX_VALUE);
		
		return new TrackGroup(annotation);
	}

	public static TrackGroup getCytobandTrackGroup(GBrowserPlot plot, CytobandDataSource cytobandData) {
		
		CytobandTrack overviewCytobands = new CytobandTrack(plot.getOverviewView(), cytobandData, false);
		return new TrackGroup(overviewCytobands);
	}

	public static TitleTrack getTitleTrack(GBrowserPlot genomePlot, String title) {
		GBrowserView dataView = genomePlot.getDataView();
		return new TitleTrack(dataView, title, Color.black);
	}
	
}
