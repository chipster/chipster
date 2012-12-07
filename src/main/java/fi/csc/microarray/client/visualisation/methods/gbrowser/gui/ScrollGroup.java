package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CoverageAndSNPTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.QualityCoverageTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackContext;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;

/**
 * ScrollGroup is a component in genome browser component hierarchy between GBrowserView and TrackGroup.
 * ScrollGroup implements vertical scrolling and therefore the size of the ScrollGroup (aka view port)
 * can be smaller than the size of its content (aka canvas). Scrolling is controlled by Swing scroll bar
 * drawn in GBrowserChartPanel. The scrolling is disabled by default.
 * 
 * Drawables are collected here from Tracks and drawn using GBrowserView's drawing implementations.
 * 
 * @author klemela
 */
public class ScrollGroup implements LayoutComponent, LayoutContainer {

	public List<TrackGroup> trackGroups = new LinkedList<TrackGroup>();
	private int layoutHeight;
	private String name;
	private boolean scrollEnabled = false;
	private BufferedImage drawBuffer;;

	/**
	 * 
	 * @param name Used only for more informative debug messages.
	 */
	public ScrollGroup(String name) {
		this.name = name;
	}

	public ScrollGroup() {
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
	 * Collects drawables from Tracks, makes last adjustments to the layout, draws the drawables and blits the right part of 
	 * the draw buffer according to scroll position. All drawables are collected before drawing any of them, because the hight of
	 * LayoutMode.FULL tracks are known only ofter the drawables are created.
	 * 
	 * @param g
	 * @param plotArea
	 * @param scrollGroupViewPort
	 * @param view
	 */
	public void draw(Graphics2D g, Rectangle plotArea, Rectangle scrollGroupViewPort, GBrowserView view) {

		/* List of all visible tracks and list of drawable collections produced by those tracks.
		 * Indexes of both lists go hand-in-hand. 
		 */
		List<Track> visibleTracks = new LinkedList<Track>();
		List<Collection<Drawable>> drawableLists = new LinkedList<Collection<Drawable>>();

		// Get drawables from all TrackGroups and Tracks
		for (TrackGroup group : getTrackGroups()) {

			Iterator<Track> trackIter = group.getTracks().iterator();
			Iterator<Drawable> drawableIter = null;

			if (!group.isVisible()) {
				continue;
			}

			// get drawables of all tracks
			while (trackIter.hasNext()) {

				Track track;

				if (drawableIter == null || !drawableIter.hasNext()) {
					track = trackIter.next();

					if (track.isVisible()) {

						Collection<Drawable> drawables = track.getDrawables();
						
						//Update track height
						track.setFullHeight(drawables);
						
						//Store drawables and Track reference
						drawableLists.add(drawables);
						visibleTracks.add(track);
					}
				}
			}
		}

		//Now we know the height of the content and can create a draw buffer for this scroll group. 
		//Ignore scrolling now and just draw everything at first.

		int drawBufferWidth = (int) (scrollGroupViewPort.getWidth());
		int drawBufferHeight = (int) (getFullHeight());

		if (drawBuffer == null || 
				drawBuffer.getWidth() != drawBufferWidth || 
				drawBuffer.getHeight() != drawBufferHeight) {		

			drawBuffer = new BufferedImage(drawBufferWidth, (int) drawBufferHeight, BufferedImage.TYPE_INT_ARGB);
		}

		Graphics2D bufG2 = (Graphics2D) drawBuffer.getGraphics();
		bufG2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);

		bufG2.setPaint(Color.white);
		bufG2.fillRect(0, 0, drawBuffer.getWidth(), drawBuffer.getHeight());

		// prepare coordinates
		int y = 0;
		int x = 0;

		Iterator<Collection<Drawable>> drawableListIter = drawableLists.iterator();
		Iterator<Track> trackIter = visibleTracks.iterator();

		// Iterate lists of track and drawables simultaneously
		while (drawableListIter.hasNext() && trackIter.hasNext()) {
			Collection<Drawable> drawables = drawableListIter.next();
			Track track = trackIter.next();
			
			//Add track height before the track drawables are drawn, because the track coordinates start
			//from the bottom and grow upwards
			y += track.getFullHeight();				

			for (Drawable drawable : drawables) {

				if(drawable == null) {
					continue;
				}

				// decide if we will expand drawable for this track
				boolean expandDrawables = track.canExpandDrawables();

				TrackContext trackContext = null;
				// create view context for this track only if we will use it
				// currently only used for tracks that contain information
				// about reads
				if (expandDrawables && 
						(track instanceof CoverageAndSNPTrack ||
								track instanceof QualityCoverageTrack)) {

					if (view.parentPlot.getReadScale() == ReadScale.AUTO) {
						trackContext = new TrackContext(track);
					} else {
						// FIXME ReadScale is in "number of reads" and context takes "number of pixels"
						trackContext = new TrackContext(track, track.getHeight() - view.parentPlot.getReadScale().numReads);
					}
				}

				// expand drawables to stretch across all height if necessary
				if (expandDrawables) {
					drawable.expand(trackContext);
				}

				// recalculate position for reversed strands
				int maybeReversedY = (int) y;
				if (track.isReversed()) {
					maybeReversedY -= track.getFullHeight();
				} else {
					drawable.upsideDown();
				}			

				// draw the drawable to the buffer
				view.drawDrawable(bufG2, x, maybeReversedY, drawable);
			}
		}               

		//Finally, get the scroll position to know which part of the content is shown
		int scrollValue = view.parentPlot.chartPanel.getScrollValue(this);
	
		int width = scrollGroupViewPort.width;
		int viewPortHeight = scrollGroupViewPort.height;

		//Fill the view port area of this ScrollGroup with the right part of the content.
		g.drawImage(drawBuffer, 
				(int) scrollGroupViewPort.getX(), 
				(int) scrollGroupViewPort.getY(), 
				(int) (scrollGroupViewPort.getX() + width), 
				(int) (scrollGroupViewPort.getY() + viewPortHeight),
				0, 
				scrollValue, 
				width, 
				scrollValue + viewPortHeight, null);

	}

	@Override
	public int getHeight() {
		return LayoutTool.getHeight(this, layoutHeight);
	}

	@Override
	public void setHeight(int height) {
		this.layoutHeight = height;
	}

	/* 
	 * Size of the content doesn't limit the minimum size of the scroll group view port,
	 * if scrolling is enabled.
	 */
	@Override
	public int getMinHeight() {
		if (scrollEnabled) {
			return 0;
		} else {
			return LayoutTool.getMinHeightSum(this);
		}
	}

	@Override
	public boolean isVisible() {
		return true;
	}

	/**
	 * Shortcut for adding a new track into this ScrollGroup. Usually related tracks should be inside the same
	 * TrackGroup, whereas this method creates always a new TrackGroup for each track. 
	 * 
	 * @param track
	 */
	public void addTrack(Track track) {
		trackGroups.add(new TrackGroup(track));
	}

	public void addTrackGroup(TrackGroup group) {
		trackGroups.add(group);
	}

	public Collection<TrackGroup> getTrackGroups() {
		return trackGroups;
	}
	
	@Override
	public Collection<? extends LayoutComponent> getLayoutComponents() {
		return getTrackGroups();
	}

	public String toString() {
		return ScrollGroup.class + " " + name;
	}

	public boolean isScrollEnabled() {
		return scrollEnabled;
	}

	@Override
	public int getFullHeight() {
		return LayoutTool.getFullHeight(this);
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

	public String getName() {
		return name;
	}

	/* 
	 * ScrollGroup view port fills the available space if scrolling is enabled.
	 */
	@Override
	public LayoutMode getLayoutMode() {
		if (scrollEnabled) {
			return LayoutMode.FILL;
		} else {
			return LayoutTool.inferLayoutMode(this);
		}
	}

	@Override
	public void setLayoutMode(LayoutMode mode) {
		// decided in getLayoutMode()
	}

	@Override
	public void setDefaultLayoutMode() {
		// decided in getLayoutMode()
	}
}
