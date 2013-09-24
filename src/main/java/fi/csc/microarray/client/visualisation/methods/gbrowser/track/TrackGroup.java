package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.OverlayLayout;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.ScrollGroup;

/**
 * <p>A collection of tracks representing a single data source or
 * related in any other logical way and displayed one after
 * another. Typically one group object corresponds to what
 * user considers as a "track".</p>
 * 
 * <p>The track group is also responsible for drawing the side menu,
 * switching tracks inside the group and changing properties of
 * individual groups.</p>
 * 
 * @author Rimvydas Naktinis, Petri Klemelä
 *
 */
public class TrackGroup implements ActionListener {
	
	JPanel component = new JPanel();
    
    protected List<Track> tracks = new LinkedList<Track>();//current tracks    
    protected GBrowserView view;

	JPanel settingsPanel;
	JCheckBox visibleBox = new JCheckBox();
	JCheckBox showMoreBox = new JCheckBox();
	private boolean isSettingsEnabled = false;
	
	private JPanel trackLayer;

	private JPanel settingsLayer;

	private Track singleTrack;

	private String name;

	private boolean minimized;

	private boolean showMore;

	private StatusAnimation statusAnimation;

	private boolean settingsPanelInitialized;

	private ScrollGroup scrollGroup;

    public TrackGroup(GBrowserView view) {
        this.view = view;
        component.setLayout(new GridBagLayout());    
        component.setInheritsPopupMenu(true);
        //Usually track fills the whole track group and listen for these events. 
        //When the height of the settings panel exceeds the height of the tracks, also the 
        //TrackGroup becomes visible and has to listen for events.
        component.addMouseListener(getView());
        component.addMouseMotionListener(getView());
        component.addMouseWheelListener(getView());
        
        component.setBackground(Color.white);     
        
		component.setLayout(new OverlayLayout(component));				
		
		trackLayer = new JPanel();
		settingsLayer = new JPanel();
		
		trackLayer.setInheritsPopupMenu(true);
		settingsLayer.setInheritsPopupMenu(true);
		
		//First component added is the topmost when drawn
		component.add(settingsLayer);
		component.add(trackLayer); 		
		
		trackLayer.setLayout(new MigLayout("flowy, fillx, gap 0! 0!, insets 0"));
		settingsLayer.setLayout(new MigLayout("gap 0! 0!, insets 0"));
		
		settingsLayer.setOpaque(false);		
		trackLayer.setBackground(Color.white);			
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
        this.singleTrack = track;
        addTrack(track);
    }
    
    /**
     * Add a track.
     */
    public void addTrack(Track track) {
        tracks.add(track);
        track.setTrackGroup(this);
        updateLayout();
    }
    
    /**
     * Return a list of tracks in this group including separator tracks.
     * 
     * @return a list of tracks do be drawn.
     */
    public List<Track> getTracks() {
        return tracks;
    }
    
    public GBrowserView getView() {
        return view;
    }      
    
    public String getName() {
    	if (name != null) {
    		return name;
    	} else {
    		return "Track Group";
    	}
    }
    
	public LayoutMode getLayoutMode() {
		return LayoutTool.inferTrackGroupLayoutMode(tracks);
	}

	public boolean isSettingsEnabled() {
		return isSettingsEnabled;
	}
	
	public void setSettingsEnabled(boolean enabled) {
		isSettingsEnabled = enabled;
	}

	private JComponent initSettingsIfNecesasry() {
		if (!settingsPanelInitialized) {
			createSettingsPanel();		
			updateButtons();
			settingsPanelInitialized = true;
		}
		return settingsPanel;
	}

	public void createSettingsPanel() {
		settingsLayer.setOpaque(false);
		
		visibleBox.setText(getName());
		Font font = visibleBox.getFont().deriveFont(Font.BOLD);
		visibleBox.setFont(font);
		visibleBox.setSelected(!isMinimized());
		visibleBox.setOpaque(false);
		visibleBox.addActionListener(this);
				
		showMoreBox.setText(getShowMoreName());
		font = showMoreBox.getFont().deriveFont(Font.PLAIN);
		showMoreBox.setFont(font);
		showMoreBox.setOpaque(false);
		showMoreBox.addActionListener(this);
				
		settingsLayer.add(visibleBox);		
		settingsLayer.add(getStatusAnimation(), "wrap");
		settingsLayer.add(showMoreBox);		
	}

	public String getShowMoreName() {
		return "Show all";
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == visibleBox) {
			setMinimzed(!visibleBox.isSelected());
		}		
		if (e.getSource() == showMoreBox) {			
			showMore(showMoreBox.isSelected());
		}
	}
	
	private void updateButtons() {
		
		boolean showMore = isShowMorePossible();
		showMoreBox.setVisible(showMore);
	}

	public boolean isShowMorePossible() {
		boolean showMore = false;
		for(Track track : tracks) {
			if (track.isVisible() && track.isSuitableViewLength() && track.isShowMoreCapable()) {
				showMore = true;
				break;
			}
		}
		return showMore;
	}
	
	protected void setMinimzed(boolean minimized) {
		this.minimized = minimized;
		addTracks();
		updateButtons();
		view.reloadData();
	}

	public void addTracks() {
		//Hook for TrackGroups to update tracks after a press of increase or decrease button
		
		tracks.clear();
		if (!isMinimized()) {
			addTrack(singleTrack);
		}
	}

	protected void showMore(boolean showMore) {
		this.showMore = showMore;
		
		addTracks();		
		updateButtons();
		view.reloadData();
	}

	public void updateLayout() {
		trackLayer.removeAll();
		
		if (isSettingsEnabled()) {
			initSettingsIfNecesasry();
		} else {
			settingsLayer.removeAll();
			settingsPanelInitialized = false;
		}
		
		for (Track track : tracks) {
			
			if (track.isVisible() && track.isSuitableViewLength()) {	        	

				Component trackComponent = track.getComponent();	        
				LayoutMode mode = track.getLayoutMode();

				if (LayoutMode.FILL == mode) {
					trackLayer.add(trackComponent, "grow");
				} else {
					trackLayer.add(trackComponent, "growx");
				}	        	
			}				
		}
		//component sizes may have changed
		trackLayer.revalidate();
		
		updateButtons();
	}

	public JComponent getComponent() {
		return component;
	}
	
	public boolean isMinimized() {
		return minimized;
	}
	
	public boolean isShowMore() {
		return showMore;
	}

	public void setName(String name) {
		this.name = name;
	}

	public StatusAnimation getStatusAnimation() {
		if (statusAnimation == null) {
			statusAnimation = new StatusAnimation(view.getQueueManager());
		}
		return statusAnimation;
	}

	public void initializeListener() {
		getStatusAnimation().initilizeListeners();
		for (Track track : tracks) {
			track.initializeListener();
		}
	}

	public ScrollGroup getScrollGroup() {
		return scrollGroup;
	}

	public void setScrollGroup(ScrollGroup scrollGroup) {
		this.scrollGroup = scrollGroup;
	}
}
