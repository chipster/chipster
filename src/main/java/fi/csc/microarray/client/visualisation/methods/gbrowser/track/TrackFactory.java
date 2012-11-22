package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowserPlot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.CytobandDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.view.View;

/**
 * Utility class for creating predefined {@link TrackGroup} objects.  
 *
 *@author Petri Klemelä, Aleksi Kallio
 *
 */
public class TrackFactory {
	
    public static TrackGroup addGeneTracks(GBrowserPlot genomePlot, DataSource annotationDataSource, TabixDataSource repeatDataSource) throws FileNotFoundException {
        
		View dataView = genomePlot.getDataView();
		
		TrackGroup geneGroup = new GeneTrackGroup(dataView, annotationDataSource, repeatDataSource);
		
		// Add gene group to data view
	    addGroup(dataView, geneGroup);
	    
	    return geneGroup;
	}

	public static void addThickSeparatorTrack(GBrowserPlot genomePlot) {
		View dataView = genomePlot.getDataView();
		dataView.addTrack(new SeparatorTrack3D(dataView, 0, Long.MAX_VALUE, false));
		dataView.addTrack(new EmptyTrack(dataView, 2));
		dataView.addTrack(new SeparatorTrack3D(dataView, 0, Long.MAX_VALUE, true));
	}
	
	public static TrackGroup addReadTracks(GBrowserPlot genomePlot, DataSource userData, DataSource seqFile, String title)
	        throws FileNotFoundException, MalformedURLException {
	
		View dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadTrackGroup readGroup = new ReadTrackGroup(dataView, userData, seqFile, title);
		readGroup.initialise();
        
        addGroup(dataView, readGroup);
        
        return readGroup;
	}
	
	public static TrackGroup addReadSummaryTracks(GBrowserPlot genomePlot, DataSource userData,
	        DataSource seqFile, String title, TabixDataSource summaryDataSource)
	        throws FileNotFoundException, MalformedURLException {
	
		View dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadSummaryTrackGroup readGroup = new ReadSummaryTrackGroup(dataView, userData, seqFile, title, summaryDataSource);
		readGroup.initialise();
		
        addGroup(dataView, readGroup);
        
        return readGroup;
	}

	public static void addWigTrack(GBrowserPlot plot, DataSource peakFile) {
		WIGTrack annotation = new WIGTrack(plot.getDataView(), peakFile,
		        Color.BLUE, 0, Long.MAX_VALUE);
		addTrack(plot.getDataView(), annotation);
	}
	
	public static void addPeakTrack(GBrowserPlot plot, DataSource peaks) {
		View dataView = plot.getDataView();

		PeakTrack annotation = new PeakTrack(dataView, peaks, GBrowserConstants.BED_COLOR, 0, Long.MAX_VALUE);
		addTrack(dataView, annotation);
	}

	public static void addTranscriptTracks(GBrowserPlot plot, DataSource transcriptAnnotationFile) {

	}

	public static void addCytobandTracks(GBrowserPlot plot, CytobandDataSource cytobandData) {
		
		CytobandTrack overviewCytobands = new CytobandTrack(plot.getOverviewView(), cytobandData, false);
		addTrack(plot.getOverviewView(), overviewCytobands);
	}

	public static void addRulerTrack(GBrowserPlot plot) {
		plot.getDataView().addTrack(new RulerTrack(plot.getDataView()));
	}
		
	private static void addTrack(View view, Track track) {
		view.addTrack(track);
		track.initializeListener();
	}
	
    private static void addGroup(View view, TrackGroup group) {
        view.addTrackGroup(group);
        
        for (Track track : group.getTracks()) {
            if (track.hasData()) {
                track.initializeListener();
            }
        }
    }

	public static void addTitleTrack(GBrowserPlot genomePlot, String title) {
		View dataView = genomePlot.getDataView();
		dataView.addTrack(new TitleTrack(dataView, title, Color.black));
	}
	
}
