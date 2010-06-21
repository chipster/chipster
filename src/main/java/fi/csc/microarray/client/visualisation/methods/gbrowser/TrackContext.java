package fi.csc.microarray.client.visualisation.methods.gbrowser;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;

/**
 * Context for a Drawable object. Context is influenced by other
 * drawable objects that are currently visible and by size of
 * the track.
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
            this.minDrawableY = Math.min(drawable.getMinY(), this.minDrawableY);
        }
               
        this.trackHeight = track.getHeight();
        
        this.expansionRatio = this.trackHeight /
                (this.trackHeight - (float) this.minDrawableY);
    }
}
