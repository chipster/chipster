package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CytobandTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.EmptyTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.GeneTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.PeakTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadSummaryTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.RulerTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TitleTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.WIGTrack;
import fi.csc.microarray.constants.VisualConstants;

/**
 * Utility class for creating predefined {@link TrackGroup} objects.  
 *
 *@author Petri Klemelä, Aleksi Kallio
 *
 */
public class TrackFactory {
	
    public static TrackGroup addGeneTracks(GenomePlot genomePlot, DataSource annotationDataSource) throws FileNotFoundException {
        
		View dataView = genomePlot.getDataView();
		
		TrackGroup geneGroup = new GeneTrackGroup(dataView, annotationDataSource);
		
		// Add gene group to data view
	    addGroup(dataView, geneGroup);
	    
	    return geneGroup;
	}

	public static void addThickSeparatorTrack(GenomePlot genomePlot) {
		View dataView = genomePlot.getDataView();
		dataView.addTrack(new SeparatorTrack3D(dataView, 0, Long.MAX_VALUE, false));
		dataView.addTrack(new EmptyTrack(dataView, 2));
		dataView.addTrack(new SeparatorTrack3D(dataView, 0, Long.MAX_VALUE, true));
	}
	
	public static TrackGroup addReadTracks(GenomePlot genomePlot, DataSource userData, DataSource seqFile, String title)
	        throws FileNotFoundException, MalformedURLException {
	
		View dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadTrackGroup readGroup = new ReadTrackGroup(dataView, userData, seqFile, title);
		readGroup.initialise();
        
        addGroup(dataView, readGroup);
        
        return readGroup;
	}
	
	public static TrackGroup addReadSummaryTracks(GenomePlot genomePlot, DataSource userData,
	        DataSource seqFile, String title, TabixDataSource summaryDataSource)
	        throws FileNotFoundException, MalformedURLException {
	
		View dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		ReadSummaryTrackGroup readGroup = new ReadSummaryTrackGroup(dataView, userData, seqFile, title, summaryDataSource);
		readGroup.initialise();
		
        addGroup(dataView, readGroup);
        
        return readGroup;
	}

	public static void addWigTrack(GenomePlot plot, DataSource peakFile) {
		WIGTrack annotation = new WIGTrack(plot.getDataView(), peakFile,
		        Color.BLUE, 0, Long.MAX_VALUE);
		addTrack(plot.getDataView(), annotation);
	}
	
	public static void addPeakTrack(GenomePlot plot, DataSource peaks) {
		View dataView = plot.getDataView();

		PeakTrack annotation = new PeakTrack(dataView, peaks, VisualConstants.BED_COLOR, 0, Long.MAX_VALUE);
		addTrack(dataView, annotation);
	}

	public static void addHeaderPeakTrack(GenomePlot plot, DataSource peaks) {
		View dataView = plot.getDataView();

		PeakTrack annotation = new PeakTrack(dataView, peaks, VisualConstants.BED_COLOR, 0, Long.MAX_VALUE);
		addTrack(dataView, annotation);
	}

	public static void addTranscriptTracks(GenomePlot plot, DataSource transcriptAnnotationFile) {

	}

	public static void addCytobandTracks(GenomePlot plot, CytobandDataSource cytobandData) {
		
		CytobandTrack overviewCytobands = new CytobandTrack(plot.getOverviewView(), cytobandData, false);
		addTrack(plot.getOverviewView(), overviewCytobands);
	}

	public static void addRulerTrack(GenomePlot plot) {
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

	public static void addTitleTrack(GenomePlot genomePlot, String title) {
		View dataView = genomePlot.getDataView();
		dataView.addTrack(new TitleTrack(dataView, title, Color.black));
	}
	
}
