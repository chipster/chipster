package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

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
public class AnnotationTrackGroup extends TrackGroup {
	
	protected TranscriptTrack transcript;
	protected Track geneOverview;
	protected Track gene;
	protected RepeatMaskerTrack repeatMasker;
	protected Track geneOverviewReversed;
	protected Track geneReversed;
	protected TranscriptTrack transcriptReversed;
	private boolean isUserData;
	private DataThread annotationDataSource;
	private DataThread repeatDataSource;
	private boolean repeat;

	public AnnotationTrackGroup(GBrowserView dataView, DataThread annotationDataSource, DataThread repeatDataSource, boolean isUserData) {
		super(dataView);
		this.isUserData = isUserData;
		this.annotationDataSource = annotationDataSource;
		this.repeatDataSource = repeatDataSource;
		initTracks();
		
		if (isUserData) {
			this.setName(annotationDataSource.getDataSource().getDataUrl().getName());
		} else {
			this.setName("Annotations");
		}
	}
	
	private void initTracks() {					
		
		getStatusAnimation().clear();
		
		if (annotationDataSource != null) {
			
			transcript = new TranscriptTrack();
			transcript.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
			transcript.setView(view);
			transcript.addDataThread(annotationDataSource);
			transcript.setStrand(Strand.FORWARD);		

			gene = new GeneTrack(GBrowserConstants.COLOR_BLUE_BRIGHTER);
			gene.setViewLimits(GBrowserConstants.SWITCH_VIEWS_AT, GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
			gene.setView(view);
			gene.addDataThread(annotationDataSource);
			gene.setStrand(Strand.FORWARD);
			
			//titleTrack.addDataThread(annotationDataSource);
			super.getStatusAnimation().addDataThread(annotationDataSource);
		}
		
		if (repeatDataSource != null) {
			
			repeatMasker = new RepeatMaskerTrack();
			repeatMasker.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
			repeatMasker.setView(view);
			repeatMasker.addDataThread(repeatDataSource);

			super.getStatusAnimation().addDataThread(repeatDataSource);
		}
		
		if (annotationDataSource != null) {

			geneReversed = new GeneTrack(GBrowserConstants.COLOR_BLUE_BRIGHTER);
			geneReversed.setViewLimits(GBrowserConstants.SWITCH_VIEWS_AT, GBrowserConstants.CHANGE_TRACKS_ZOOM_THRESHOLD2);
			geneReversed.setView(view);
			geneReversed.addDataThread(annotationDataSource);
			geneReversed.setStrand(Strand.REVERSE);

			transcriptReversed = new TranscriptTrack();
			transcriptReversed.setViewLimits(0, GBrowserConstants.SWITCH_VIEWS_AT);
			transcriptReversed.setView(view);
			transcriptReversed.addDataThread(annotationDataSource);
			transcriptReversed.setStrand(Strand.REVERSE);
		}
		
		addTracks();
	}

	@Override
	public void addTracks() {
		
		tracks.clear();
						
        if (!isMinimized() && transcript != null) { // no annotation data source 
        	// Transcript, detailed, forward
        	addTrack(transcript);

        	// Gene, detailed, forward
        	addTrack(gene);
        }

		if (isUserData) {
			if (!isMinimized()) {
				SeparatorTrack separator = new SeparatorTrack(Color.gray, 1);
				separator.setView(view);
				addTrack(separator);
			}
		} else {
			
			if (isMinimized()) {
				SeparatorTrack separator = new SeparatorTrack(Color.white, 20);
				separator.setView(view);
				addTrack(separator);
			}
			
			// Ruler track
			RulerTrack ruler = new RulerTrack();
			ruler.setView(view);
			addTrack(ruler);			
		}

		if (!isMinimized() && repeatMasker != null) {
			if (repeat) {
				// Repeat masker track
				addTrack(repeatMasker);
			}
		}
		
		if (!isMinimized() && transcript != null) { //no annotation data source

			// Gene, detailed, reverse
			addTrack(geneReversed);

			// Transcript, detailed, reverse
			addTrack(transcriptReversed);
		}
		
		if (isShowMore()) {
			if (gene != null && geneReversed != null) {
				gene.setLayoutMode(LayoutMode.FULL);
				geneReversed.setLayoutMode(LayoutMode.FULL);
			}

			if (transcript != null && transcriptReversed != null) {
				transcript.setLayoutMode(LayoutMode.FULL);
				transcriptReversed.setLayoutMode(LayoutMode.FULL);
			}
		} else {
			if (gene != null && geneReversed != null) {
				gene.setDefaultLayoutMode();
				geneReversed.setDefaultLayoutMode();
			}
			if (transcript != null && transcriptReversed != null) {
				transcript.setDefaultLayoutMode();
				transcriptReversed.setDefaultLayoutMode();
			}
		}
	}

	public void setRepeatVisible(boolean selected) {
		this.repeat = selected;
		addTracks();
	}
}
