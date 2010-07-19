package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;

import fi.csc.microarray.client.visualisation.methods.gbrowser.View;

/**
 * A collection of tracks representing a single data source or
 * related in any other logical way and displayed one after
 * another.
 *  
 * @author naktinis
 *
 */
public class TrackGroup {
    
    protected List<Track> tracks = new LinkedList<Track>();
    protected View view;
    protected boolean menuVisible;
    public JPanel menu;
    
    // Width of side menu in pixels
    public static final int MENU_WIDTH = 70;

    public TrackGroup(View view) {
        this.view = view;
        
        // Add side menu
        menu = new JPanel();
        this.view.parentPlot.chartPanel.add(menu);
    }
    
    /**
     * It is quite common to have track group consisting of a
     * single data track.
     * 
     * @param view
     * @param track
     */
    public TrackGroup(Track track) {
        this(track.view);
        tracks.add(track);
    }
    
    /**
     * Add a track.
     */
    public void addTrack(Track track) {
        tracks.add(track);
    }
    
    /**
     * Return a list of tracks in this group including separator tracks.
     * 
     * @return a list of tracks do be drawn.
     */
    public List<Track> getTracks() {
        // TODO add separator tracks
        return tracks;
    }
    
    public int getHeight() {
        int height = 0;
        for (Track track : getTracks()) {
            height += track.getHeight();
        }
        return height;
    }
    
    public View getView() {
        return view;
    }
    
    /**
     * Determine if a side menu should be shown for this track.
     * 
     * @return true if a menu should be shown, false otherwise.
     */
    public boolean isMenuVisible() {
        return menuVisible;
    }
    
    /**
     * Set side menu visibility. 
     */
    public void setMenuVisible(boolean isVisible) {
        menuVisible = isVisible;
    }
}
