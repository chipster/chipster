package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CytobandTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.GeneTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.IntensityTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.PeakTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.RepeatMaskerTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.RulerTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TitleTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.WIGTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack.PartColor;

/**
 * Utility for creating predefined track groups.  
 *
 */
public class TrackFactory {
	
	private static final int CHANGE_TRACKS_ZOOM_THRESHOLD2 = 10000000;
	private static final int CHANGE_TRACKS_ZOOM_THRESHOLD1 = 100000;

    public static TrackGroup addGeneTracks(GenomePlot genomePlot, ChunkDataSource geneAnnotationFile,
	        DataSource transcriptAnnotationFile, ChunkDataSource refSource) throws FileNotFoundException {
        
		View dataView = genomePlot.getDataView();
		
		// Group containing gene-related tracks
        TrackGroup geneGroup = new TrackGroup(dataView);
        
        // Top separator and title
        geneGroup.addTrack(createThickSeparatorTrack(dataView));
        geneGroup.addTrack(new TitleTrack(dataView, "Annotations", Color.black));
		
		// Transcript, detailed, forward
		TranscriptTrack trancsript = new TranscriptTrack(dataView, transcriptAnnotationFile, ChunkTreeHandlerThread.class,
		        Color.DARK_GRAY, CHANGE_TRACKS_ZOOM_THRESHOLD1);
		trancsript.setStrand(Strand.FORWARD);
		geneGroup.addTrack(trancsript);

		// Gene, overview, forward 
		IntensityTrack annotationOverview = new IntensityTrack(genomePlot.getDataView(),
		        geneAnnotationFile, ChunkTreeHandlerThread.class, PartColor.CDS.c, CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotationOverview.setStrand(Strand.FORWARD);
		geneGroup.addTrack(annotationOverview);

		// Gene, detailed, forward
		GeneTrack annotation = new GeneTrack(genomePlot.getDataView(), geneAnnotationFile,
		        ChunkTreeHandlerThread.class, PartColor.CDS.c, CHANGE_TRACKS_ZOOM_THRESHOLD1, CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotation.setStrand(Strand.FORWARD);
		geneGroup.addTrack(annotation);
			
		// Ruler track
		geneGroup.addTrack(new RulerTrack(dataView));
		  
        // Repeat masker track
        RepeatMaskerTrack repeatMasker =
            new RepeatMaskerTrack(dataView, refSource, ChunkTreeHandlerThread.class, Long.MAX_VALUE);
        geneGroup.addTrack(repeatMasker);
		
		// Gene, overview, reverse
		IntensityTrack annotationOverviewReversed = new IntensityTrack(genomePlot.getDataView(),
		        geneAnnotationFile, ChunkTreeHandlerThread.class, PartColor.CDS.c, CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotationOverviewReversed.setStrand(Strand.REVERSED);
        geneGroup.addTrack(annotationOverviewReversed);

		// Gene, detailed, reverse
		GeneTrack annotationReversed = new GeneTrack(genomePlot.getDataView(), geneAnnotationFile,
		        ChunkTreeHandlerThread.class, PartColor.CDS.c, CHANGE_TRACKS_ZOOM_THRESHOLD1, CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotationReversed.setStrand(Strand.REVERSED);
        geneGroup.addTrack(annotationReversed);
	
		// Transcript, detailed, reverse
		TranscriptTrack trancsriptReverse = new TranscriptTrack(dataView, transcriptAnnotationFile, ChunkTreeHandlerThread.class,
		        Color.DARK_GRAY, CHANGE_TRACKS_ZOOM_THRESHOLD1);
		trancsriptReverse.setStrand(Strand.REVERSED);
		geneGroup.addTrack(trancsriptReverse);
		
		// Add gene group to data view
	    addGroup(dataView, geneGroup);
	    
	    return geneGroup;
	}

	private static void addSeparatorTrack(GenomePlot genomePlot) {
		addSeparatorTrack(genomePlot, Long.MAX_VALUE);
	}
	
	private static void addSeparatorTrack(GenomePlot genomePlot, long maxBpLength) {
		View dataView = genomePlot.getDataView();
		dataView.addTrack(new SeparatorTrack(dataView, Color.gray, 1, 0, maxBpLength));
	}

	static void addThickSeparatorTrack(GenomePlot genomePlot) {
		View dataView = genomePlot.getDataView();
		dataView.addTrack(createThickSeparatorTrack(dataView));
	}
	
    public static Track createThickSeparatorTrack(View view) {
        return new SeparatorTrack(view, Color.gray.brighter(), 4, 0, Long.MAX_VALUE);
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

	// FIXME Currently not used, used miRNAParser
	public static void addMirnaTracks(GenomePlot genomePlot, ChunkDataSource miRNAFile) {
		View dataView = genomePlot.getDataView();

		for (Strand strand : Strand.values()) {

			GeneTrack track = new GeneTrack(dataView, miRNAFile, ChunkTreeHandlerThread.class,
			        PartColor.CDS.c.darker(), 0, Long.MAX_VALUE);
			track.setStrand(strand);
			dataView.addTrack(track);
			track.initializeListener();

			if (strand == Strand.FORWARD) {
				addSeparatorTrack(genomePlot);
			}
		}
	}

	public static void addCytobandTracks(GenomePlot plot, ChunkDataSource cytobandData) {
		CytobandTrack overviewCytobands = new CytobandTrack(plot.getOverviewView(), cytobandData, ChunkTreeHandlerThread.class, false);
		addTrack(plot.getOverviewView(), overviewCytobands);

		//CytobandTrack cytobands = new CytobandTrack(plot.getDataView(), cytobandData, ChunkTreeHandlerThread.class, true);
		//addTrack(plot.getDataView(), cytobands);
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
