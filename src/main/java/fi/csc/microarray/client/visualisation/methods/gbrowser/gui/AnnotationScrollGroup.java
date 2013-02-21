package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.GeneTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TranscriptTrack;


public class AnnotationScrollGroup extends ScrollGroup {

	public AnnotationScrollGroup() {
		super("Annotations", true);
	}
	
	@Override
	public int getScrollReferenceY() {
		//Return the height of the first GeneTrack or TranscriptTrack, i.e. keep the 
		//RulerTrack steady when the height of the genes or transcripts changes.
		for (Track track : trackGroups.get(0).getTracks()) {
			if (track instanceof GeneTrack || track instanceof TranscriptTrack) {
				return Math.max(track.getFullHeight(), track.getHeight()); 
			}
		}
		return 0;
	}
	
	@Override
	public LayoutMode getLayoutMode() {
		return LayoutTool.inferLayoutMode(this);
	}
}
