package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

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
	private StatusTitleTrack titleTrack;

	public GeneTrackGroup(GBrowserView dataView, DataThread annotationDataSource, DataThread repeatDataSource, boolean isUserData) {
		super(dataView);
		
		titleTrack = new StatusTitleTrack("Annotations", Color.black);
		titleTrack.setView(view);
		
		if (annotationDataSource != null) {
			
			transcript = new TranscriptTrack();
			transcript.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
			transcript.setView(dataView);
			transcript.addDataThread(annotationDataSource);
			transcript.setStrand(Strand.FORWARD);
			

//			geneOverview = new CoverageEstimateTrack(dataView, annotationDataSource, GBrowserConstants.COLOR_BLUE_BRIGHTER, 
//					GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true, false);
//			geneOverview.setStrand(Strand.FORWARD);
//			geneOverview = new EmptyTrack(transcript.getMinHeight(), GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
//			geneOverview.setView(dataView);

			gene = new GeneTrack(GBrowserConstants.COLOR_BLUE_BRIGHTER);
			gene.setViewLimits(GBrowserConstants.SWITCH_VIEWS_AT, GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
			gene.setView(dataView);
			gene.addDataThread(annotationDataSource);
			gene.setStrand(Strand.FORWARD);
			
			titleTrack.addDataThread(annotationDataSource);
		}
		
		if (repeatDataSource != null) {
			repeatMasker = new RepeatMaskerTrack();
			repeatMasker.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
			repeatMasker.setView(dataView);
			repeatMasker.addDataThread(repeatDataSource);
			
			titleTrack.addDataThread(repeatDataSource);
		}
		
		if (annotationDataSource != null) {
//			geneOverviewReversed = new CoverageEstimateTrack(dataView, annotationDataSource, GBrowserConstants.COLOR_BLUE_BRIGHTER, 
//					GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2, true, false);
//			geneOverviewReversed.setStrand(Strand.REVERSE);
//			geneOverviewReversed = new EmptyTrack(transcript.getMinHeight(), GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
//			geneOverviewReversed.setView(dataView);

			geneReversed = new GeneTrack(GBrowserConstants.COLOR_BLUE_BRIGHTER);
			geneReversed.setViewLimits(GBrowserConstants.SWITCH_VIEWS_AT, GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
			geneReversed.setView(dataView);
			geneReversed.addDataThread(annotationDataSource);
			geneReversed.setStrand(Strand.REVERSE);

			transcriptReversed = new TranscriptTrack();
			transcriptReversed.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
			transcriptReversed.setView(dataView);
			transcriptReversed.addDataThread(annotationDataSource);
			transcriptReversed.setStrand(Strand.REVERSE);
		}
		
		adds(isUserData);
	}

	public void adds(boolean isUserData) {
		
		this.tracks = new LinkedList<Track>();
		
		if (!isUserData) {
			// title
			tracks.add(titleTrack);
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
			SeparatorTrack separator = new SeparatorTrack(Color.gray, 1);
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
