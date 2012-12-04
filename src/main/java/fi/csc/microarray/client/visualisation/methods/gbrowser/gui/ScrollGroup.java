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
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CoverageAndSNPTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.QualityCoverageTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackContext;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;

public class ScrollGroup implements LayoutComponent, LayoutContainer {

	public List<TrackGroup> trackGroups = new LinkedList<TrackGroup>();
	private int layoutHeight;
	private ScrollPosition defaultScrollPosition = ScrollPosition.MID;
	private String name;
	//	private int maxViewPortHeight = 300;
	private boolean scrollEnabled = false;
	private BufferedImage drawBuffer;;

	public ScrollGroup(String name) {
		this.name = name;
	}

	public ScrollGroup(String name, ScrollPosition defaultScrollPosition) {
		this(name);
		this.defaultScrollPosition = defaultScrollPosition;
	}

	//	public ScrollGroup(String name, ScrollPosition defaultScrollPosition, int maxViewPortHeight) {
	//		this(name, defaultScrollPosition);
	//		this.maxViewPortHeight = maxViewPortHeight;
	//	}

	public ScrollGroup() {
	}

	public ScrollGroup(String name, boolean scrollEnabled) {
		this(name);
		this.scrollEnabled = scrollEnabled;
	}

	public void draw(Graphics2D g, Rectangle plotArea, Rectangle scrollGroupViewPort, GBrowserView view) {

		// prepare context object
		List<Collection<Drawable>> drawableLists = new LinkedList<Collection<Drawable>>();
		List<Track> visibleTracks = new LinkedList<Track>();

		// track group contains one or several logically-related tracks
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

					// draw drawable objects for visible tracks
					if (track.isVisible()) {

						Collection<Drawable> drawables = track.getDrawables();
						track.setFullHeight(drawables);
						drawableLists.add(drawables);
						visibleTracks.add(track);
					}
				}
			}
		}


		int drawBufferWidth = (int) (scrollGroupViewPort.getWidth());
		int drawBufferHeight = (int) (getFullHeight());

		if (drawBuffer == null || 
				drawBuffer.getWidth() != drawBufferWidth || 
				drawBuffer.getHeight() != drawBufferHeight) {		

			/* drawBuffer contains only this view. Plot coordinates have to be shifted by the size of the other views
			 * (viewArea.x and viewArea.y)
			 */
			drawBuffer = new BufferedImage(drawBufferWidth, (int) drawBufferHeight, BufferedImage.TYPE_INT_ARGB);
		}

		Graphics2D bufG2 = (Graphics2D) drawBuffer.getGraphics();
		bufG2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
				RenderingHints.VALUE_ANTIALIAS_ON);

		/* The JScrollPane doesn't clip the content properly, but it is drawn outside JScrollPane when the window is resized.
		 * Probably this has something to do with our custom use of clip areas or maybe the JFreeChart uses some ancient AWT 
		 * components. Nevertheless, making the content transparent and drawing only the JViewPort area solves the problem, 
		 * as the drawing transparent pixels elsewhere doesn't have any effect. 
		 */
		//			bufG2.setBackground(new Color(0, 0, 0, 0));			

		//		if (isFullHeight && !isFixedHeight()) {
		//			
		//			//bufG2.setClip(null);
		//			bufG2.clearRect(0, 0, drawBuffer.getWidth(), drawBuffer.getHeight());
		//
		//			
		//			/* In full height mode we try to draw always also the content that isn't shown in current vertical
		//			 * scrolling position. Setting the original clip to the drawing buffer should
		//			 * at least prevent actual pixel manipulating when drawing outside of the view.
		//			 */
		//			Rectangle clipRectangle = new Rectangle(plotArea.x - viewArea.x, plotArea.y - viewArea.y, plotArea.width, plotArea.height);
		//			bufG2.setClip(clipRectangle);
		//						
		//			bufG2.setPaint(Color.white);
		//			bufG2.fill(clipRectangle);
		//			
		//		} else {
		//			bufG2.setClip(null);

		bufG2.setPaint(Color.white);
		bufG2.fillRect(0, 0, drawBuffer.getWidth(), drawBuffer.getHeight());
		//		}

		// prepare coordinates
		int y = 0;
		int x = 0;

		Iterator<Collection<Drawable>> drawableListIter = drawableLists.iterator();
		Iterator<Track> trackIter = visibleTracks.iterator();

		while (drawableListIter.hasNext() && trackIter.hasNext()) {

			Collection<Drawable> drawables = drawableListIter.next();
			Track track = trackIter.next();
			
			//Add track height before the track drawables are drawn, because the track coordinates start
			//from the bottom and grow upwards
			y += track.getFullHeight();				


			//			int maxY = 20;

			//						if (isFullHeight && track.isFixedHeight()) {
			//
			//							while (drawableIter.hasNext()) {										
			//
			//								Drawable drawable = drawableIter.next();
			//
			//								if(drawable == null) {
			//									continue;
			//								}
			//
			//								if (drawable.getMaxY() > maxY) {
			//									maxY = drawable.getMaxY();
			//								}
			//							}
			//
			//							track.setCanvasHeight(maxY + FULL_HEIGHT_MARGIN);						
			//						}


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

				// draw an object on the buffer
				view.drawDrawable(bufG2, x, maybeReversedY, drawable);
			}
		}               


		int scrollValue = view.parentPlot.chartPanel.getScrollValue(this);

		//copy only the visible area of the JScrollPane
		//			g.drawImage(drawBuffer, 
		//					(int) plotViewPort.getX(), 
		//					(int) plotViewPort.getY() + scrollGroupDestinationY, 
		//					(int) (plotViewPort.getX() + plotViewPort.getWidth()), 
		//					(int) (plotViewPort.getY() + scrollGroup.getHeight() + scrollGroupDestinationY),
		//					(int) plotViewPort.getX() - viewCanvasArea.x, 
		//					(int) plotViewPort.getY() - viewCanvasArea.y + scrollGroupCanvasY + scrollValue, 
		//					(int) (plotViewPort.getX() - viewCanvasArea.x + plotViewPort.getWidth()), 
		//					(int) (plotViewPort.getY() - viewCanvasArea.y + scrollGroup.getHeight() + scrollGroupCanvasY  + scrollValue), null);
		
		int width = scrollGroupViewPort.width;
		int viewPortHeight = scrollGroupViewPort.height;

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
	public Collection<? extends LayoutComponent> getLayoutComponents() {
		return trackGroups;
	}

	@Override
	public int getHeight() {
		return LayoutTool.getHeight(this, layoutHeight);
	}

	@Override
	public void setHeight(int height) {
		this.layoutHeight = height;
	}

	@Override
	public int getMinHeight() {
		return 0;
	}

	@Override
	public boolean isVisible() {
		return true;
	}

	public void addTrack(Track track) {
		trackGroups.add(new TrackGroup(track));
	}

	public void addTrackGroup(TrackGroup group) {
		trackGroups.add(group);
	}

	public Collection<TrackGroup> getTrackGroups() {
		return trackGroups;
	}

	public ScrollPosition getDefaultScrollPosition() {
		return defaultScrollPosition ;
	}

	public enum ScrollPosition { START, MID };

	public String toString() {
		return ScrollGroup.class + " " + name;
	}

	//FIXME can be removed?
	//	public int getMaxViewPortHeight() {
	//		return maxViewPortHeight ;
	//	}

	public boolean isScrollEnabled() {
		return scrollEnabled;
	}

	@Override
	public int getFullHeight() {
		return LayoutTool.getFullHeight(this);
	}

	public int getScrollReferenceY() {
		return 0;
	}

	public String getName() {
		return name;
	}

	@Override
	public LayoutMode getLayoutMode() {
		if (scrollEnabled) {
			return LayoutMode.FILL;
		} else {
			return LayoutTool.inferLayoutMode(this);
		}
	}

	public void setFullLayoutMode(boolean enabled) {
		for (TrackGroup group : trackGroups) {
			group.setFullLayoutMode(enabled);
		}
	}
}
