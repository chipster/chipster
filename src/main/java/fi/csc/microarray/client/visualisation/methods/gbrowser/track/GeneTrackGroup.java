package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.constants.VisualConstants;

/**
 * Track group containing information about genes: transcript, intensity, gene, snp
 * repeat masker.
 * 
 * @author Vilius Zukauskas
 *
 */
public class GeneTrackGroup extends TrackGroup {
	
	protected TranscriptTrack transcript;
	protected IntensityTrack geneOverview;
	protected Track gene;
	protected ReferenceSNPTrack snpTrack = null;
	protected RepeatMaskerTrack repeatMasker;
	protected IntensityTrack geneOverviewReversed;
	protected Track geneReversed;
	protected TranscriptTrack transcriptReversed;
	protected ReferenceSNPTrack snpTrackReversed;

	public GeneTrackGroup(View dataView, ChunkDataSource geneAnnotationFile,
	        DataSource transcriptAnnotationFile, ChunkDataSource refSource, DataSource snpFile) {
		super(dataView);
		
		transcript = new TranscriptTrack(dataView, transcriptAnnotationFile, ChunkTreeHandlerThread.class,
		        Color.DARK_GRAY, GenomeBrowserConstants.SWITCH_VIEWS_AT);
		transcript.setStrand(Strand.FORWARD);
		
		geneOverview = new IntensityTrack(dataView, geneAnnotationFile, ChunkTreeHandlerThread.class, 
				VisualConstants.COLOR_BLUE_BRIGHTER, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true, false);
		geneOverview.setStrand(Strand.FORWARD);
		
		gene = new GeneTrack(dataView, geneAnnotationFile,
		        ChunkTreeHandlerThread.class, VisualConstants.COLOR_BLUE_BRIGHTER, GenomeBrowserConstants.SWITCH_VIEWS_AT, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
		gene.setStrand(Strand.FORWARD);
		
		if (snpFile != null) {
			snpTrack = new ReferenceSNPTrack(dataView, snpFile, ChunkTreeHandlerThread.class, 0, GenomeBrowserConstants.SHOW_SNP_AT);
			snpTrack.setStrand(Strand.FORWARD);

			snpTrackReversed = new ReferenceSNPTrack(dataView, snpFile, ChunkTreeHandlerThread.class, 0, GenomeBrowserConstants.SHOW_SNP_AT);
			snpTrackReversed.setStrand(Strand.REVERSED);
		}
		
		repeatMasker = new RepeatMaskerTrack(dataView, refSource, ChunkTreeHandlerThread.class, GenomeBrowserConstants.SWITCH_VIEWS_AT);
		
		geneOverviewReversed = new IntensityTrack(dataView,
		        geneAnnotationFile, ChunkTreeHandlerThread.class, VisualConstants.COLOR_BLUE_BRIGHTER, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true, false);
		geneOverviewReversed.setStrand(Strand.REVERSED);
		
		geneReversed = new GeneTrack(dataView, geneAnnotationFile,
		        ChunkTreeHandlerThread.class, VisualConstants.COLOR_BLUE_BRIGHTER, GenomeBrowserConstants.SWITCH_VIEWS_AT, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
		geneReversed.setStrand(Strand.REVERSED);
		
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
		tracks.add(geneOverview);

		// Gene, detailed, forward
		tracks.add(gene);
		
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
        tracks.add(geneOverviewReversed);

		// Gene, detailed, reverse
        tracks.add(geneReversed);
	
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
