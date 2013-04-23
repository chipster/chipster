package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;

/**
 * Track group containing information about genes: transcript, intensity, gene, snp
 * repeat masker.
 * 
 * @author Vilius Zukauskas, Petri Klemelä
 *
 */
public class GeneTrackGroup extends TrackGroup {
	
	protected TranscriptTrack transcript;
	protected Track geneOverview;
	protected Track gene;
	protected RepeatMaskerTrack repeatMasker;
	protected Track geneOverviewReversed;
	protected Track geneReversed;
	protected TranscriptTrack transcriptReversed;

	public GeneTrackGroup(GBrowserView dataView, AreaRequestHandler annotationDataSource, AreaRequestHandler repeatDataSource, boolean isUserData) {
		super(dataView);
		
		if (annotationDataSource != null) {
			transcript = new TranscriptTrack(GBrowserConstants.SWITCH_VIEWS_AT);
			transcript.setView(dataView);
			transcript.addAreaRequestHandler(annotationDataSource);
			transcript.setStrand(Strand.FORWARD);

//			geneOverview = new CoverageEstimateTrack(dataView, annotationDataSource, GBrowserConstants.COLOR_BLUE_BRIGHTER, 
//					GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true, false);
//			geneOverview.setStrand(Strand.FORWARD);
//			geneOverview = new EmptyTrack(transcript.getMinHeight(), GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
//			geneOverview.setView(dataView);

			gene = new GeneTrack(GBrowserConstants.COLOR_BLUE_BRIGHTER, 
					GBrowserConstants.SWITCH_VIEWS_AT, GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
			gene.setView(dataView);
			gene.addAreaRequestHandler(annotationDataSource);
			gene.setStrand(Strand.FORWARD);
		}
		
		if (repeatDataSource != null) {
			repeatMasker = new RepeatMaskerTrack(0, GBrowserConstants.SWITCH_VIEWS_AT);
			repeatMasker.setView(dataView);
			repeatMasker.addAreaRequestHandler(repeatDataSource);
		}
		
		if (annotationDataSource != null) {
//			geneOverviewReversed = new CoverageEstimateTrack(dataView, annotationDataSource, GBrowserConstants.COLOR_BLUE_BRIGHTER, 
//					GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true, false);
//			geneOverviewReversed.setStrand(Strand.REVERSE);
//			geneOverviewReversed = new EmptyTrack(transcript.getMinHeight(), GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
//			geneOverviewReversed.setView(dataView);

			geneReversed = new GeneTrack(GBrowserConstants.COLOR_BLUE_BRIGHTER, 
					GBrowserConstants.SWITCH_VIEWS_AT, GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
			geneReversed.setView(dataView);
			geneReversed.addAreaRequestHandler(annotationDataSource);
			geneReversed.setStrand(Strand.REVERSE);

			transcriptReversed = new TranscriptTrack(GBrowserConstants.SWITCH_VIEWS_AT);
			transcriptReversed.setView(dataView);
			transcriptReversed.addAreaRequestHandler(annotationDataSource);
			transcriptReversed.setStrand(Strand.REVERSE);
		}
		
		adds(isUserData);
	}

	public void adds(boolean isUserData) {
		
		this.tracks = new LinkedList<Track>();
		
		if (!isUserData) {
			// title
			
			TitleTrack title = new TitleTrack("Annotations", Color.black);
			title.setView(view);
			tracks.add(title);
		}
		
        if (transcript != null) { // no annotation data source 
        	// Transcript, detailed, forward
        	tracks.add(transcript);

        	// Gene, overview, forward 
        	//tracks.add(geneOverview);

        	// Gene, detailed, forward
        	tracks.add(gene);
        }

		if (isUserData) {
			SeparatorTrack separator = new SeparatorTrack(Color.gray, 1, 0, Long.MAX_VALUE);
			separator.setView(view);
			tracks.add(separator);
		} else {
			// Ruler track
			RulerTrack ruler = new RulerTrack();
			ruler.setView(view);
			tracks.add(ruler);			
		}

		if (repeatMasker != null) {
			// Repeat masker track
			tracks.add(repeatMasker);
		}
		
		if (transcript != null) { //no annotation data source
			// Gene, overview, reverse
			//tracks.add(geneOverviewReversed);

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
	
	@Override
	public void showOrHide(String name, boolean state) {
		super.showOrHide(name, state);
	}
	
	@Override
	public LayoutMode getLayoutMode() {
		return LayoutMode.FIXED;
	}
	
	@Override
	public int getHeight() {
		if ( getView().getBpRegion().getLength() <  GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2) {
			return 250;
		} else {
			return 40;
		}
	}	
}
