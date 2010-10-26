package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CytobandTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.GeneTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.PeakTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadSummaryTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.RulerTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TitleTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.WIGTrack;

/**
 * Utility for creating predefined track groups.  
 *
 */
public class TrackFactory {
	
    public static TrackGroup addGeneTracks(GenomePlot genomePlot, ChunkDataSource geneAnnotationFile,
	        DataSource transcriptAnnotationFile, ChunkDataSource refSource, DataSource snpFile) throws FileNotFoundException {
        
		View dataView = genomePlot.getDataView();
		
		TrackGroup geneGroup = new GeneTrackGroup(dataView, geneAnnotationFile,
				transcriptAnnotationFile, refSource, snpFile);
		
		// Add gene group to data view
	    addGroup(dataView, geneGroup);
	    
	    return geneGroup;
	}

	static void addThickSeparatorTrack(GenomePlot genomePlot) {
		View dataView = genomePlot.getDataView();
		dataView.addTrack(createThickSeparatorTrack(dataView));
	}
	
    public static Track createThickSeparatorTrack(View view) {
        return new SeparatorTrack(view, Color.gray.brighter(), true, 0, Long.MAX_VALUE);
    }

	public static TrackGroup addReadTracks(GenomePlot genomePlot, DataSource userData,
	        Class<? extends AreaRequestHandler> userDataHandler,
	        DataSource seqFile, String title)
	        throws FileNotFoundException, MalformedURLException {
	
		View dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		TrackGroup readGroup = new ReadTrackGroup(dataView, userData,
		        userDataHandler, seqFile, title);
        
        addGroup(dataView, readGroup);
        
        return readGroup;
	}
	
	public static TrackGroup addReadSummaryTracks(GenomePlot genomePlot, DataSource userData,
	        Class<? extends AreaRequestHandler> userDataHandler, DataSource seqFile, String filename)
	        throws FileNotFoundException, MalformedURLException {
	
		View dataView = genomePlot.getDataView();
		
		// Group containing tracks for this data source
		TrackGroup readGroup = new ReadSummaryTrackGroup(dataView, userData,
		        userDataHandler, seqFile, filename);
        
        addGroup(dataView, readGroup);
        
        return readGroup;
	}

	public static void addWigTrack(GenomePlot plot, DataSource peakFile) {
		WIGTrack annotation = new WIGTrack(plot.getDataView(), peakFile,
		        ChunkTreeHandlerThread.class, Color.BLUE, 0, Long.MAX_VALUE);
		addTrack(plot.getDataView(), annotation);
	}
	
	public static void addPeakTrack(GenomePlot plot, DataSource peaks) {
		View dataView = plot.getDataView();

		PeakTrack annotation = new PeakTrack(dataView, peaks, ChunkTreeHandlerThread.class, Color.YELLOW, 0, Long.MAX_VALUE);
		addTrack(dataView, annotation);
	}

	public static void addHeaderPeakTrack(GenomePlot plot, DataSource peaks) {
		View dataView = plot.getDataView();

		PeakTrack annotation = new PeakTrack(dataView, peaks, ChunkTreeHandlerThread.class, Color.YELLOW, 0, Long.MAX_VALUE);
		addTrack(dataView, annotation);
	}

	public static void addTranscriptTracks(GenomePlot plot, DataSource transcriptAnnotationFile) {

	}

	public static void addCytobandTracks(GenomePlot plot, ChunkDataSource cytobandData) {
		CytobandTrack overviewCytobands = new CytobandTrack(plot.getOverviewView(), cytobandData, ChunkTreeHandlerThread.class, false);
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
