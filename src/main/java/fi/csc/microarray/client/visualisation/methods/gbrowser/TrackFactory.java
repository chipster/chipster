package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CytobandTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.GelTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.GeneTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.IntensityTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.PeakTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ProfileTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.RulerTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeqBlockTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeqTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TitleTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack.PartColor;

/**
 * Utility for creating predefined track groups.  
 *
 */
public class TrackFactory {
	
	public static void addGeneTracks(GenomePlot genomePlot, ChunkDataSource geneAnnotationFile,
	        DataSource transcriptAnnotationFile) {

		View dataView = genomePlot.getDataView();
		
		// Gene, overview, forward 
		IntensityTrack annotationOverview = new IntensityTrack(genomePlot.getDataView(),
		        geneAnnotationFile, ChunkTreeHandlerThread.class, PartColor.CDS.c, 10000000);
		annotationOverview.setStrand(Strand.FORWARD);
		addTrack(dataView, annotationOverview);

		// Gene, detailed, forward
		GeneTrack annotation = new GeneTrack(genomePlot.getDataView(), geneAnnotationFile,
		        ChunkTreeHandlerThread.class, PartColor.CDS.c, 0, 10000000);
		annotation.setStrand(Strand.FORWARD);
		addTrack(dataView, annotation);
		
		// Add Transcript tracks for both strands
		for (Strand strand : Strand.values()) {

			// Transcript, overview
			IntensityTrack transcriptOverview = new IntensityTrack(dataView, transcriptAnnotationFile,
			        ChunkTreeHandlerThread.class, PartColor.CDS.c.darker(), 100000);
			transcriptOverview.setStrand(strand);
			addTrack(dataView, transcriptOverview);

			// Transcript, detailed
			TranscriptTrack trancsript = new TranscriptTrack(dataView, transcriptAnnotationFile,
			        ChunkTreeHandlerThread.class, Color.DARK_GRAY, 100000);
			trancsript.setStrand(strand);
			addTrack(dataView, trancsript);

			if (strand == Strand.FORWARD) {
				addSeparatorTrack(genomePlot);
			}
		}
		
		// Gene, overview, reverse 
		IntensityTrack annotationOverviewReversed = new IntensityTrack(genomePlot.getDataView(),
		        geneAnnotationFile, ChunkTreeHandlerThread.class, PartColor.CDS.c, 10000000);
		annotationOverviewReversed.setStrand(Strand.REVERSED);
		addTrack(dataView, annotationOverviewReversed);

		// Gene, detailed, reverse
		GeneTrack annotationReversed = new GeneTrack(genomePlot.getDataView(), geneAnnotationFile,
		        ChunkTreeHandlerThread.class, PartColor.CDS.c, 0, 10000000);
		annotationReversed.setStrand(Strand.REVERSED);
		addTrack(dataView, annotationReversed);
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
		dataView.addTrack(new SeparatorTrack(dataView, Color.gray.brighter(), 4, 0, Long.MAX_VALUE));
	}

	public static void addReadTracks(GenomePlot genomePlot, DataSource userData,
	        Class<? extends AreaRequestHandler> userDataHandler,
	        DataSource seqFile, boolean isTreatment)
	        throws FileNotFoundException, MalformedURLException {
	
		View dataView = genomePlot.getDataView();
		int switchViewsAt = 50000;
		Color histogramColor = isTreatment ? Color.blue : Color.gray;
		Color fontColor = Color.black;
							
		// 
		// Forward
		//

		// Overview
		IntensityTrack readOverview = new IntensityTrack(dataView, userData, userDataHandler, histogramColor, switchViewsAt);
		addTrack(dataView, readOverview);

		// Detailed
		SeqBlockTrack reads = new SeqBlockTrack(dataView, userData, userDataHandler, fontColor, 0, switchViewsAt);
		addTrack(dataView, reads);
		
	    // Gel
        addSeparatorTrack(genomePlot);
        GelTrack reads2 = new GelTrack(dataView, userData, userDataHandler, Color.GREEN, 0, Long.MAX_VALUE);
        addTrack(dataView, reads2);

		// Profile
	    addSeparatorTrack(genomePlot);
		ProfileTrack reads3 = new ProfileTrack(dataView, userData, userDataHandler, fontColor, 0, Long.MAX_VALUE);
		addTrack(dataView, reads3);

		addSeparatorTrack(genomePlot);

		//
		// Reference sequence
		//

		if (seqFile != null) {
			// Reference sequence
			SeqTrack seq = new SeqTrack(dataView, seqFile, ChunkTreeHandlerThread.class, 800);
			addTrack(dataView, seq);
			addSeparatorTrack(genomePlot, 800);
		}

		//
		// Reverse
		//

		// Overview
		IntensityTrack readOverviewReversed = new IntensityTrack(dataView, userData, userDataHandler, histogramColor, switchViewsAt);
		readOverviewReversed.setStrand(Strand.REVERSED);
		addTrack(dataView, readOverviewReversed);

		// Detailed
		SeqBlockTrack readsReversed = new SeqBlockTrack(dataView, userData, userDataHandler, fontColor, 0, switchViewsAt);
		readsReversed.setStrand(Strand.REVERSED);
		addTrack(dataView, readsReversed);
	}

	// FIXME Currently not used, used miRNAParser
	public static void addWigTrack(GenomePlot plot, DataSource peakFile) {
		ProfileTrack annotation = new ProfileTrack(plot.getDataView(), peakFile,
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

		CytobandTrack cytobands = new CytobandTrack(plot.getDataView(), cytobandData, ChunkTreeHandlerThread.class, true);
		addTrack(plot.getDataView(), cytobands);
	}

	public static void addRulerTrack(GenomePlot plot) {
		plot.getDataView().addTrack(new RulerTrack(plot.getDataView()));
	}
	
	private static void addTrack(View view, Track track) {
		view.addTrack(track);
		track.initializeListener();
	}

	public static void addTitleTrack(GenomePlot genomePlot, String title) {
		View dataView = genomePlot.getDataView();
		dataView.addTrack(new TitleTrack(dataView, title, Color.black));
	}
	

}
