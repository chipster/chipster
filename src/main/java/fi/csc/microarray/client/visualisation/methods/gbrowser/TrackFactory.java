package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TSVHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.CytobandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.GeneParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.HeaderTsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SequenceParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TranscriptParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.miRNAParser;
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
	
	private static final int CHANGE_TRACKS_ZOOM_THRESHOLD2 = 10000000;
	private static final int CHANGE_TRACKS_ZOOM_THRESHOLD1 = 100000;

	public static void addGeneTracks(GenomePlot genomePlot, DataSource geneAnnotationFile, DataSource transcriptAnnotationFile) {

		// initialise data source files
		GeneParser geneParser = new GeneParser();
		TranscriptParser transcriptParser = new TranscriptParser();
		View dataView = genomePlot.getDataView();
		
		// Transcript, detailed, forward
		TranscriptTrack trancsript = new TranscriptTrack(dataView, transcriptAnnotationFile, TSVHandlerThread.class, transcriptParser, Color.DARK_GRAY, CHANGE_TRACKS_ZOOM_THRESHOLD1);
		trancsript.setStrand(Strand.FORWARD);
		addTrack(dataView, trancsript);

		// Gene, overview, forward 
		IntensityTrack annotationOverview = new IntensityTrack(genomePlot.getDataView(), geneAnnotationFile, TSVHandlerThread.class, geneParser, PartColor.CDS.c, CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotationOverview.setStrand(Strand.FORWARD);
		addTrack(dataView, annotationOverview);

		// Gene, detailed, forward
		GeneTrack annotation = new GeneTrack(genomePlot.getDataView(), geneAnnotationFile, TSVHandlerThread.class, geneParser, PartColor.CDS.c, CHANGE_TRACKS_ZOOM_THRESHOLD1, CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotation.setStrand(Strand.FORWARD);
		addTrack(dataView, annotation);
		
		// Separator line
		addRulerTrack(genomePlot);

		// Gene, overview, reverse 
		IntensityTrack annotationOverviewReversed = new IntensityTrack(genomePlot.getDataView(), geneAnnotationFile, TSVHandlerThread.class, geneParser, PartColor.CDS.c, CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotationOverviewReversed.setStrand(Strand.REVERSED);
		addTrack(dataView, annotationOverviewReversed);

		// Gene, detailed, reverse
		GeneTrack annotationReversed = new GeneTrack(genomePlot.getDataView(), geneAnnotationFile, TSVHandlerThread.class, geneParser, PartColor.CDS.c, CHANGE_TRACKS_ZOOM_THRESHOLD1, CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotationReversed.setStrand(Strand.REVERSED);
		addTrack(dataView, annotationReversed);

		// Transcript, detailed, reverse
		TranscriptTrack trancsriptReverse = new TranscriptTrack(dataView, transcriptAnnotationFile, TSVHandlerThread.class, transcriptParser, Color.DARK_GRAY, CHANGE_TRACKS_ZOOM_THRESHOLD1);
		trancsriptReverse.setStrand(Strand.REVERSED);
		addTrack(dataView, trancsriptReverse);
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


	public static void addReadTracks(GenomePlot genomePlot, DataSource userData, DataSource seqFile, boolean isTreatment) throws FileNotFoundException, MalformedURLException {
		addReadTracks(genomePlot, userData, seqFile, isTreatment, new ElandParser());
	}

	public static void addReadTracks(GenomePlot genomePlot, DataSource userData, DataSource seqFile, boolean isTreatment, TsvParser userDataParser) throws FileNotFoundException, MalformedURLException {
	
		View dataView = genomePlot.getDataView();
		int switchViewsAt = 50000;
		Color histogramColor = isTreatment ? Color.blue : Color.gray;
		Color fontColor = Color.black;
							
		// 
		// Forward
		//

		// Overview
		IntensityTrack readOverview = new IntensityTrack(dataView, userData, TSVHandlerThread.class, userDataParser, histogramColor, switchViewsAt);
		addTrack(dataView, readOverview);

		// Detailed
		SeqBlockTrack reads = new SeqBlockTrack(dataView, userData, TSVHandlerThread.class, userDataParser, fontColor, 0, switchViewsAt);
		addTrack(dataView, reads);
		
	    // Gel
        addSeparatorTrack(genomePlot);
        GelTrack reads2 = new GelTrack(dataView, userData, TSVHandlerThread.class, userDataParser, Color.GREEN, 0, Long.MAX_VALUE);
        addTrack(dataView, reads2);

		// Profile
	    addSeparatorTrack(genomePlot);
		ProfileTrack reads3 = new ProfileTrack(dataView, userData, TSVHandlerThread.class, userDataParser, fontColor, 0, Long.MAX_VALUE);
		addTrack(dataView, reads3);

		addSeparatorTrack(genomePlot);

		//
		// Reference sequence
		//

		if (seqFile != null) {
			// Reference sequence
			SeqTrack seq = new SeqTrack(dataView, seqFile, TSVHandlerThread.class, new SequenceParser(), 800);
			addTrack(dataView, seq);
			addSeparatorTrack(genomePlot, 800);
		}

		//
		// Reverse
		//

		// Overview
		IntensityTrack readOverviewReversed = new IntensityTrack(dataView, userData, TSVHandlerThread.class, userDataParser, histogramColor, switchViewsAt);
		readOverviewReversed.setStrand(Strand.REVERSED);
		addTrack(dataView, readOverviewReversed);

		// Detailed
		SeqBlockTrack readsReversed = new SeqBlockTrack(dataView, userData, TSVHandlerThread.class, userDataParser, fontColor, 0, switchViewsAt);
		readsReversed.setStrand(Strand.REVERSED);
		addTrack(dataView, readsReversed);
	}

	public static void addWigTrack(GenomePlot plot, DataSource peakFile) {
		miRNAParser miRNAParser = new miRNAParser();
		ProfileTrack annotation = new ProfileTrack(plot.getDataView(), peakFile, TSVHandlerThread.class, miRNAParser, Color.BLUE, 0, Long.MAX_VALUE);
		addTrack(plot.getDataView(), annotation);
	}
	
	public static void addPeakTrack(GenomePlot plot, DataSource peaks) {
		BEDParser bedParser = new BEDParser();
		View dataView = plot.getDataView();

		PeakTrack annotation = new PeakTrack(dataView, peaks, TSVHandlerThread.class, bedParser, Color.YELLOW, 0, Long.MAX_VALUE);
		addTrack(dataView, annotation);
	}

	public static void addHeaderPeakTrack(GenomePlot plot, DataSource peaks) {
		HeaderTsvParser headerTsvParser = new HeaderTsvParser();
		View dataView = plot.getDataView();

		PeakTrack annotation = new PeakTrack(dataView, peaks, TSVHandlerThread.class, headerTsvParser, Color.YELLOW, 0, Long.MAX_VALUE);
		addTrack(dataView, annotation);
	}

	public static void addTranscriptTracks(GenomePlot plot, DataSource transcriptAnnotationFile) {

	}

	public static void addMirnaTracks(GenomePlot genomePlot, DataSource miRNAFile) {
		miRNAParser miRNAParser = new miRNAParser();
		View dataView = genomePlot.getDataView();

		for (Strand strand : Strand.values()) {

			GeneTrack track = new GeneTrack(dataView, miRNAFile, TSVHandlerThread.class, miRNAParser, PartColor.CDS.c.darker(), 0, Long.MAX_VALUE);
			track.setStrand(strand);
			dataView.addTrack(track);
			track.initializeListener();

			if (strand == Strand.FORWARD) {
				addSeparatorTrack(genomePlot);
			}
		}
	}

	public static void addCytobandTracks(GenomePlot plot, DataSource cytobandFile) {
		CytobandTrack overviewCytobands = new CytobandTrack(plot.getOverviewView(), cytobandFile, TSVHandlerThread.class, new CytobandParser(), false);
		addTrack(plot.getOverviewView(), overviewCytobands);

		CytobandTrack cytobands = new CytobandTrack(plot.getDataView(), cytobandFile, TSVHandlerThread.class, new CytobandParser(), true);
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
