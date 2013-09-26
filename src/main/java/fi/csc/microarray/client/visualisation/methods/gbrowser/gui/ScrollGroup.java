package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JPanel;
import javax.swing.JScrollPane;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;

/**
 * ScrollGroup is a component in genome browser component hierarchy between GBrowserView and TrackGroup.
 * ScrollGroup implements vertical scrolling and therefore the size of the ScrollGroup (aka view port)
 * can be smaller than the size of its content (aka canvas). The scrolling is disabled by default.
 * 
 * Drawables are collected here from Tracks and drawn using GBrowserView's drawing implementations.
 * 
 * @author klemela
 */
public class ScrollGroup {
	
	private JScrollPane component = new JScrollPane() {
		@Override
		public Dimension getMinimumSize() {
			//Set minimum height to preferred height for small tracks, because there is not point to scroll them
			//Minimum height for larger scroll groups is 250 pixels, because smaller view port is not that useful 
			//either
			Dimension pref = super.getPreferredSize();			
			
			if (pref.getHeight() < 250) {
				return pref;
			}
			return new Dimension((int) pref.getWidth(), 250);
			
		}
				
		@Override
		public Dimension getMaximumSize() {
			//Do not let view port to grow larger than canvas, unless the canvas is able to fill that space			
			Dimension pref = canvas.getPreferredSize();
			Dimension max = super.getMaximumSize();
			
			if (LayoutMode.FILL == getLayoutMode()) {
				return max;
			} else {
				return pref;
			}
		}
	};

	public List<TrackGroup> trackGroups = new LinkedList<TrackGroup>();
	private JPanel canvas = new JPanel();
	private String name;
	private boolean scrollEnabled = false;

	private GBrowserView view;

	public ScrollGroup() {
		super();
		component.setViewportView(canvas);
		component.getVerticalScrollBar().setUnitIncrement(10);
		component.setBorder(null);
		canvas.setLayout(new MigLayout("flowy, fillx, gap 0! 0!, insets 0"));
		component.setInheritsPopupMenu(true);
		canvas.setInheritsPopupMenu(true);
				
		/*
		 * The cursor line in GBrowserView has to be redrawn when the scroll pane is scrolled.
		 * This is not a perfect solution, because the component is shown for a while before
		 * the cursor line is painted (the line flickers during the scrolling).
		 */
		component.getVerticalScrollBar().addAdjustmentListener(new AdjustmentListener() {			
			@Override
			public void adjustmentValueChanged(AdjustmentEvent e) {
				if (view != null){					
					view.getComponent().repaint();
				}
			}
		});
	}
	
	/**
	 * 
	 * @param name Used only for more informative debug messages.
	 */
	public ScrollGroup(String name) {
		this();
		this.name = name;
	}

	/**
	 * Use this constructor to enable vertical scrolling.
	 * 
	 * @param name Used only for more informative debug messages.
	 * @param scrollEnabled
	 */
	public ScrollGroup(String name, boolean scrollEnabled) {
		this(name);
		this.scrollEnabled = scrollEnabled;
	}

	/**
	 * Shortcut for adding a new track into this ScrollGroup. Usually related tracks should be inside the same
	 * TrackGroup, whereas this method creates always a new TrackGroup for each track. 
	 * 
	 * @param track
	 */
	public void addTrack(Track track) {
		addTrackGroup(new TrackGroup(track));
	}

	public void addTrackGroup(TrackGroup group) {
		trackGroups.add(group);
		group.setScrollGroup(this);
		updateLayout();
	}

	public Collection<TrackGroup> getTrackGroups() {
		return trackGroups;
	}

	public String toString() {
		return ScrollGroup.class + " " + name;
	}

	public boolean isScrollEnabled() {
		return scrollEnabled;
	}

	/**
	 * Vertical position of the content that is kept steady when the size of the content changes.
	 * Default value is 0, which keeps the top part of the content visible.
	 * 
	 * @return
	 */
	public int getScrollReferenceY() {
		return 0;
	}	

	public String getScrollGroupName() {
		return name;
	}

	public LayoutMode getLayoutMode() {
		
		LayoutMode inferedMode = LayoutTool.inferScrollGroupLayoutMode(trackGroups);	
		return inferedMode;
	}

	public void setScrollEnabled(boolean enabled) {
		this.scrollEnabled = enabled;
	}

	public void updateLayout() {
		
		canvas.removeAll();
		
		for (TrackGroup group : trackGroups) {
			group.updateLayout();
				        
	        LayoutMode mode = group.getLayoutMode();
	        
	        if (LayoutMode.FIXED == mode) {
	        	canvas.add(group.getComponent(), "growx");
	        } else { //LayoutMode.FILL
	        	canvas.add(group.getComponent(), "grow");
	        }
		}
		
		//scroll pane canvas size may have changed
		canvas.revalidate();
	}

	public Container getCanvas() {
		return canvas;
	}

	public JScrollPane getComponent() {
		return component;
	}
	
	public GBrowserView getView() {
		return view;
	}
	
	public void setView(GBrowserView view) {
		this.view = view;
	}
}
