package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack.PartColor;

/**
 * track group containing information about genes: transcript, intensity, gene, snp
 * repeat masker.
 * 
 * @author zukauska
 *
 */

public class GeneTrackGroup extends TrackGroup {
	
	protected TranscriptTrack transcript;
	protected IntensityTrack annotationOverview;
	protected GeneTrack annotation;
	protected ReferenceSNPTrack snpTrack = null;
	protected RepeatMaskerTrack repeatMasker;
	protected IntensityTrack annotationOverviewReversed;
	protected GeneTrack annotationReversed;
	protected TranscriptTrack transcriptReversed;
	protected ReferenceSNPTrack snpTrackReversed;

	public GeneTrackGroup(View dataView, ChunkDataSource geneAnnotationFile,
	        DataSource transcriptAnnotationFile, ChunkDataSource refSource, DataSource snpFile) {
		super(dataView);
		
		transcript = new TranscriptTrack(dataView, transcriptAnnotationFile, ChunkTreeHandlerThread.class,
		        Color.DARK_GRAY, GenomeBrowserConstants.SWITCH_VIEWS_AT);
		transcript.setStrand(Strand.FORWARD);
		
		annotationOverview = new IntensityTrack(dataView, geneAnnotationFile, ChunkTreeHandlerThread.class, 
				PartColor.CDS.c, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true);
		annotationOverview.setStrand(Strand.FORWARD);
		
		annotation = new GeneTrack(dataView, geneAnnotationFile,
		        ChunkTreeHandlerThread.class, PartColor.CDS.c, GenomeBrowserConstants.SWITCH_VIEWS_AT, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotation.setStrand(Strand.FORWARD);
		
		if (snpFile != null) {
			snpTrack = new ReferenceSNPTrack(dataView, snpFile, ChunkTreeHandlerThread.class, 0, GenomeBrowserConstants.SHOW_SNP_AT);
			snpTrack.setStrand(Strand.FORWARD);

			snpTrackReversed = new ReferenceSNPTrack(dataView, snpFile, ChunkTreeHandlerThread.class, 0, GenomeBrowserConstants.SHOW_SNP_AT);
			snpTrackReversed.setStrand(Strand.REVERSED);
		}
		
		repeatMasker = new RepeatMaskerTrack(dataView, refSource, ChunkTreeHandlerThread.class, GenomeBrowserConstants.SWITCH_VIEWS_AT);
		
		annotationOverviewReversed = new IntensityTrack(dataView,
		        geneAnnotationFile, ChunkTreeHandlerThread.class, PartColor.CDS.c, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true);
		annotationOverviewReversed.setStrand(Strand.REVERSED);
		
		annotationReversed = new GeneTrack(dataView, geneAnnotationFile,
		        ChunkTreeHandlerThread.class, PartColor.CDS.c, GenomeBrowserConstants.SWITCH_VIEWS_AT, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
		annotationReversed.setStrand(Strand.REVERSED);
		
		transcriptReversed = new TranscriptTrack(dataView, transcriptAnnotationFile, ChunkTreeHandlerThread.class,
		        Color.DARK_GRAY, GenomeBrowserConstants.SWITCH_VIEWS_AT);
		transcriptReversed.setStrand(Strand.REVERSED);
		
		adds();
	}

	public void adds() {
		
		this.tracks = new LinkedList<Track>();
		// Top separator and title
        tracks.add(new TitleTrack(view, "Annotations", Color.black));
		
		// Transcript, detailed, forward
		
		tracks.add(transcript);

		// Gene, overview, forward 
		tracks.add(annotationOverview);

		// Gene, detailed, forward
		tracks.add(annotation);
		
		if (snpTrack != null) {
			// SNP track Forward
			tracks.add(snpTrack);
		}

		// Ruler track
		tracks.add(new RulerTrack(view));

		if (snpTrackReversed != null) {
			// SNP track Reversed
			tracks.add(snpTrackReversed);
		}
		  
        // Repeat masker track
        tracks.add(repeatMasker);
		
		// Gene, overview, reverse
        tracks.add(annotationOverviewReversed);

		// Gene, detailed, reverse
        tracks.add(annotationReversed);
	
		// Transcript, detailed, reverse
		tracks.add(transcriptReversed);
		
		// Add gene group to data view
//	    addGroup(view, tracks);
	}
	
	@Override
	public String getName() {
		return "GeneTrackGroup";
	}
	
	private void setChangeSNP(boolean change) {
		if (change) {
			snpTrack.changeSNPView(ChunkTreeHandlerThread.class);
			snpTrackReversed.changeSNPView(ChunkTreeHandlerThread.class);
		} else {
			snpTrack.returnSNPView();
			snpTrackReversed.returnSNPView();
		}
		view.fireAreaRequests();
        view.redraw();
	}
	
	@Override
	public void showOrHide(String name, boolean state) {
		super.showOrHide(name, state);
		if (snpTrack != null && name.equals("changeSNP")) {
			setChangeSNP(state);
		}
	}
}
