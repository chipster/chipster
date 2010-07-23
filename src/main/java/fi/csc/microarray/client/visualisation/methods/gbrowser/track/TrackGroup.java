package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedList;
import java.util.List;

import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComponent;
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
    protected boolean menuVisible = false;
    public SideMenu menu;
    private JButton resize;
        
    public class SideMenu extends JPanel implements ActionListener {
        
        // Width of side menu in pixels
        public static final int WIDTH = 110;
        public static final int COLLAPSED_WIDTH = 20;
        
        // Position of this menu
        int x, y;
        
        // Inner panel for holding controls
        private JPanel controls = new JPanel();
        int startControlsAt = 25;
        int lastControlAt = 0;
        
        protected boolean menuCollapsed = true;
        protected TrackGroup group = TrackGroup.this;
        
        public SideMenu() {
            // Absolute positioning inside the menu
            this.setLayout(null);
            try {
                resize = new JButton(
                        new ImageIcon(this.getClass().
                        getResource("/arrow_left.png").toURI().toURL()));
            } catch (Exception e) {
                e.printStackTrace();
            }
            resize.setBackground(new Color(0, 0, 0, 0));
            resize.setBorder(null);
            resize.setBounds(1, 1, 20, 20);
            resize.setOpaque(false);
            resize.setFocusPainted(false);
            resize.addActionListener(this);

            // Panel with controls
            controls.setLayout(new BoxLayout(controls, BoxLayout.PAGE_AXIS));
            controls.setBounds(5, startControlsAt, 100, 100);
            controls.setVisible(!menuCollapsed);
            add(controls);
            
            // Button for resizing side menu
            add(resize);
        }
        
        public void addItem(JComponent component) {
            controls.add(component);
        }
        
        public int getWidth() {
            return menuCollapsed ? COLLAPSED_WIDTH : WIDTH;
        }
        
        /**
         * Sets the absolute position and draws the menu.
         * 
         * @param x rightmost pixel of the menu
         * @param y topmost pixel of the menu
         */
        public void setPosition(int x, int y) {
            this.x = x;
            this.y = y;
            this.setBounds(x - getWidth(), y, getWidth(), group.getHeight());
        }
        
        /**
         * Redraw menu after collapsing.
         */
        public void redraw() {
            controls.setVisible(!menuCollapsed);
            setPosition(x, y);
        }

        @Override
        public void actionPerformed(ActionEvent event) {
            if (event.getSource() == resize) {
                menuCollapsed = !menuCollapsed;
                redraw();
            }
        }
    }

    public TrackGroup(View view) {
        this.view = view;
        
        // Add side menu
        menu = new SideMenu();
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
