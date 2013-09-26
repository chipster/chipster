package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;

/**
 * Helper class for collecting drawing related information for a {@link Drawable} 
 * object. Context is influenced by other
 * drawable objects that are currently visible and by size of
 * the track. It is used when placing and sizing {@link Drawable} objects to a {@link Track}.
 * 
 * @author naktinis
 *
 */
public class TrackContext {
    
    // Maximum y coordinate of a drawable element in this context 
    private Integer minDrawableY;
    
    // Height of the track in which visualisation takes place
    public Integer trackHeight;
    
    // Ratio at which drawables in this context can be expanded/contracted
    public Float expansionRatio;
    
    public TrackContext(Track track) {
        this.minDrawableY = Integer.MAX_VALUE;
        
        for (Drawable drawable : track.getDrawables()) {
            // Drawables can only be drawn between 0 and height-1
            this.minDrawableY = Math.min(drawable.getMaxY() + 1, this.minDrawableY);
        }
                      
        this.trackHeight = track.getTrackHeight();
        
        this.expansionRatio = this.trackHeight /
                (this.trackHeight - (float) this.minDrawableY);
    }
    
    public TrackContext(Track track, Integer minDrawableY) {           
        this.trackHeight = track.getTrackHeight();
        
        this.minDrawableY = minDrawableY;
        
//        this.expansionRatio = this.trackHeight /
//                (this.trackHeight - (float) this.minDrawableY);
        
        this.expansionRatio = this.trackHeight / (float) this.minDrawableY;
    }
}
