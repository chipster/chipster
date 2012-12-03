package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.constants.VisualConstants;

/**
 * Track group containing information about genes: transcript, intensity, gene, snp
 * repeat masker.
 * 
 * @author Vilius Zukauskas, Petri Klemelä
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

	public GeneTrackGroup(View dataView, DataSource annotationDataSource, TabixDataSource repeatDataSource) {
		super(dataView);
		
		if (annotationDataSource != null) {
			transcript = new TranscriptTrack(dataView, annotationDataSource, GenomeBrowserConstants.SWITCH_VIEWS_AT);
			transcript.setStrand(Strand.FORWARD);

			geneOverview = new IntensityTrack(dataView, annotationDataSource, VisualConstants.COLOR_BLUE_BRIGHTER, 
					GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true, false);
			geneOverview.setStrand(Strand.FORWARD);

			gene = new GeneTrack(dataView, annotationDataSource, VisualConstants.COLOR_BLUE_BRIGHTER, 
					GenomeBrowserConstants.SWITCH_VIEWS_AT, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
			gene.setStrand(Strand.FORWARD);
		}
		
		if (repeatDataSource != null) {
			repeatMasker = new RepeatMaskerTrack(dataView, repeatDataSource, 0, GenomeBrowserConstants.SWITCH_VIEWS_AT);
		}
		
		if (annotationDataSource != null) {
			geneOverviewReversed = new IntensityTrack(dataView, annotationDataSource, VisualConstants.COLOR_BLUE_BRIGHTER, 
					GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true, false);
			geneOverviewReversed.setStrand(Strand.REVERSED);

			geneReversed = new GeneTrack(dataView, annotationDataSource, VisualConstants.COLOR_BLUE_BRIGHTER, 
					GenomeBrowserConstants.SWITCH_VIEWS_AT, GenomeBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
			geneReversed.setStrand(Strand.REVERSED);

			transcriptReversed = new TranscriptTrack(dataView, annotationDataSource, GenomeBrowserConstants.SWITCH_VIEWS_AT);
			transcriptReversed.setStrand(Strand.REVERSED);
		}
		
		adds();
	}

	public void adds() {
		
		this.tracks = new LinkedList<Track>();
		// Top separator and title
        tracks.add(new TitleTrack(view, "Annotations", Color.black));
		
        if (transcript != null) { // no annotation data source 
        	// Transcript, detailed, forward
        	tracks.add(transcript);

        	// Gene, overview, forward 
        	tracks.add(geneOverview);

        	// Gene, detailed, forward
        	tracks.add(gene);
        }
		
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

		if (repeatMasker != null) {
			// Repeat masker track
			tracks.add(repeatMasker);
		}
		
		if (transcript != null) { //no annotation data source
			// Gene, overview, reverse
			tracks.add(geneOverviewReversed);

			// Gene, detailed, reverse
			tracks.add(geneReversed);

			// Transcript, detailed, reverse
			tracks.add(transcriptReversed);
		}
		
		// Add gene group to data view
//	    addGroup(view, tracks);
	}
	
	@Override
	public String getName() {
		return "GeneTrackGroup";
	}
	
	private void setChangeSNP(boolean change) {
		if (change) {
			snpTrack.changeSNPView();
			snpTrackReversed.changeSNPView();
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
