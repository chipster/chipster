package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.geom.Point2D;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.SwingUtilities;
import javax.swing.Timer;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TooltipAugmentedChartPanel.TooltipRequestProcessor;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRetrievalStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.QueueManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;

/**
 * Combines tracks to create a single synchronised view. All tracks within one view move horizontally hand-in-hand. 
 * View is responsible for handling user input that moves or zooms the view. 
 * Also the drawing of drawables is implemented in View, because different views 
 * (e.g. Horizontal or Circular) require different drawing implementations.
 * 
 */
public abstract class GBrowserView implements MouseListener, MouseMotionListener, MouseWheelListener, TooltipRequestProcessor, LayoutComponent, LayoutContainer{

	public RegionDouble bpRegion;
	public Region highlight;

	public Collection<ScrollGroup> scrollGroups = new LinkedList<ScrollGroup>();
	protected Rectangle viewArea = new Rectangle(0, 0, 500, 500);
	private QueueManager queueManager;
	private Point2D dragStartPoint;
	private boolean dragStarted;

	public GBrowserPlot parentPlot;

	protected static final int FPS = 30;

	protected boolean movable;
	protected boolean zoomable;

	protected final float ZOOM_FACTOR = 1.06f;
	protected final float MIN_PIXELS_PER_NUCLEOTIDE = 10f;

	private List<RegionListener> listeners = new LinkedList<RegionListener>();
	public int margin = 0;
	private Point2D dragEndPoint;
	private Point2D dragLastStartPoint;
	private long dragEventTime;

	private static final long DRAG_EXPIRATION_TIME_MS = 50;

	public GBrowserView(GBrowserPlot parent, boolean movable, boolean zoomable, boolean selectable) {
		this.parentPlot = parent;
		this.movable = movable;
		this.zoomable = zoomable;
	}

	/**
	 * Override to implement drawing of drawables.
	 * 
	 * @param g
	 * @param x
	 * @param y
	 * @param drawable
	 */
	protected abstract void drawDrawable(Graphics2D g, int x, int y, Drawable drawable);

	/**
	 * Add a track group containing one or several tracks.
	 * 
	 * @param group
	 */
	public void addTrackGroup(TrackGroup group) {
		ScrollGroup scrollGroup = new ScrollGroup();
		scrollGroup.addTrackGroup(group);
		scrollGroups.add(scrollGroup);
	}

	/**
	 * Get tracks contained in track groups as a linear collection.
	 * 
	 * Only return tracks within visible groups.
	 */
	public Collection<Track> getTracks() {
		Collection<Track> tracks = new LinkedList<Track>();
		for (ScrollGroup scrollGroup : scrollGroups) {
			for (TrackGroup trackGroup : scrollGroup.getTrackGroups()) {
				// Only return tracks within visible groups
				if (trackGroup.isVisible()) {
					tracks.addAll(trackGroup.getTracks());
				}
			}
		}
		return tracks;
	}

	/**
	 * Calls ScrollGroups to draw themselves.
	 * 
	 * @param g
	 * @param plotArea
	 * @param viewArea
	 */
	public void draw(Graphics2D g, Rectangle plotArea, Rectangle viewArea) {
		
		this.viewArea = (Rectangle) viewArea.clone();
				
		if (viewArea.height == 0) {
			//There is no layout yet
			return;
		}

		if (bpRegion == null) {
			setBpRegion(new RegionDouble(0d, 1024 * 1024 * 250d, new Chromosome("1")), false);
		}
		
		
		Rectangle scrollGroupViewPort = (Rectangle) viewArea.clone();
		
		for (ScrollGroup scrollGroup : scrollGroups) {
			
			scrollGroupViewPort.height = scrollGroup.getHeight();
			
			scrollGroup.draw(g, plotArea, scrollGroupViewPort, this);	
			
			scrollGroupViewPort.y += scrollGroupViewPort.height;
		}
		
		// draw simple vertical cursor line
		if (isCursorLineEnabled()) {
			g.setColor(new Color(0f, 0f, 0f, 0.25f));
			g.drawLine((int)plotArea.getWidth()/2, 0, (int)plotArea.getWidth()/2, (int)plotArea.getHeight());
		}
		
		// Print fps
//		frameCount++;
//		if (System.currentTimeMillis() >= resetTime + 500) {
//			
//			System.out.println("FPS:\t" + frameCount * 1000 / (System.currentTimeMillis() - resetTime));
//			resetTime = System.currentTimeMillis();
//			frameCount = 0;
//		}
	}
	
	// Used above
//	private int frameCount = 0;
//	private long resetTime = 0;

	public boolean isCursorLineEnabled() {
		return true;
	}

	public int getWidth() {
		return this.viewArea.width;
	}

	/**
	 * Height of the area drawn by this view. The height is sum of component heights if the components have fixed height,
	 * otherwise height set by LayoutTool.doLayout(). 
	 * 
	 * @return   
	 */
	public int getHeight() {
		return LayoutTool.getHeight(this, layoutHeight);
	}

	public QueueManager getQueueManager() {	
		if (queueManager == null) {
			queueManager = new QueueManager();
		}
		return queueManager;
	}

	/**
	 * Fire area requests for all tracks in this view.
	 * 
	 * Only fire one request for a single file.
	 */
	public void fireAreaRequests() {

		Map<AreaRequestHandler, Set<ColumnType>> datas = new HashMap<AreaRequestHandler, Set<ColumnType>>();

		// Add all requested columns for each requested file
		for (Track t : getTracks()) {
			Map<AreaRequestHandler, Set<ColumnType>> trackDatas = t.requestedData();

			// Don't do anything for hidden tracks or tracks without data
			if (trackDatas == null || !t.isVisible()) {
				continue;
			}

			for (AreaRequestHandler file : trackDatas.keySet()) {
				
				if (file != null) {
					// Add columns for this requested file
					Set<ColumnType> columns = datas.get(file);
					columns = columns != null ? columns : new HashSet<ColumnType>();
					columns.addAll(trackDatas.get(file));
					datas.put(file, columns);
				}
			}
		}
		
		Region requestRegion = getBpRegion();

		// Fire area requests
		for (AreaRequestHandler file : datas.keySet()) {
			DataRetrievalStatus status = new DataRetrievalStatus();
			status.clearQueues = true;
			getQueueManager().addAreaRequest(file, new AreaRequest(requestRegion, datas.get(file), status), true);
		}
	}
	
	public long getMinBp(long length) {
		return (long) (-length / 30);
	}
	
	public long getMaxBp(long length) {
		return (long) (length * (1 + 1.0 / 30));
	}

	public void setBpRegion(RegionDouble region, boolean disableDrawing) {				

		this.bpRegion = region;
		
		limitRegion();

		fireAreaRequests();
		
		dispatchRegionChange();
	}

	public RegionDouble getBpRegionDouble() {
		return bpRegion;
	}

	public Region getBpRegion() {
			return new Region((long) (double) bpRegion.start.bp, bpRegion.start.chr, (long)Math.ceil((double) bpRegion.end.bp), bpRegion.end.chr);
	}

	public void mouseClicked(MouseEvent e) {
		// ignore
	}

	public void mouseEntered(MouseEvent e) {
		// ignore
	}

	public void mouseExited(MouseEvent e) {
		// ignore
	}

	public void mousePressed(MouseEvent e) {

		parentPlot.chartPanel.requestFocusInWindow();

		stopAnimation();
		dragStartPoint = scale(e.getPoint());
		
		dragStarted = false;
	}

	public void mouseReleased(MouseEvent e) {

		if (dragStarted && dragEndPoint != null && dragLastStartPoint != null && Math.abs(dragEndPoint.getX() - dragLastStartPoint.getX()) > 10 && System.currentTimeMillis() - dragEventTime < DRAG_EXPIRATION_TIME_MS) {

			stopAnimation();

			mouseAnimationTimer = new Timer(1000 / FPS, new ActionListener() {

				private int i = 2; //Skip a few frames to get a head start
				private int ANIMATION_FRAMES = 30;
				private long startTime = System.currentTimeMillis();

				public void actionPerformed(ActionEvent arg0) {
					
					boolean skipFrame = false;
					boolean done = false;

					do {
						double endX = dragEndPoint.getX();
						double startX = dragLastStartPoint.getX();

						double newX = endX - (endX - startX) / (ANIMATION_FRAMES - i);

						dragEndPoint = new Point2D.Double(newX, dragEndPoint.getY());

						skipFrame = (i < (ANIMATION_FRAMES - 1)) && System.currentTimeMillis() > startTime + (1000 / FPS) * i;

						done = i >= ANIMATION_FRAMES;
						
						if (!done) {
							handleDrag(dragLastStartPoint, dragEndPoint, skipFrame);
							i++;
						} else {
							stopAnimation();
						}
						
					} while (skipFrame && !done);
				}
			});
			mouseAnimationTimer.setCoalesce(true);
			mouseAnimationTimer.setRepeats(true);
			mouseAnimationTimer.start();
		}
	}

	public void mouseDragged(MouseEvent e) {

		if (movable && ((dragStartPoint != null && viewArea.contains(dragStartPoint) || viewArea.contains(e.getPoint())))) {

			dragStarted = true;
			dragEndPoint = scale(e.getPoint());
			dragEventTime = System.currentTimeMillis();

			handleDrag(dragStartPoint, dragEndPoint, false);

		}
		dragLastStartPoint = dragStartPoint;
		dragStartPoint = scale(e.getPoint());
	}

	public void mouseMoved(MouseEvent e) {
		// ignore
	}

	protected abstract void handleDrag(Point2D start, Point2D end, boolean disableDrawing);

	private Timer mouseAnimationTimer;

	private ViewLimiter viewLimiter;
	private int layoutHeight;

	public void mouseWheelMoved(final MouseWheelEvent e) {

		zoomAnimation((int) scale(e.getPoint()).getX(), e.getWheelRotation());
	}

	public void zoomAnimation(final int centerX, final int wheelRotation) {
		stopAnimation();

		mouseAnimationTimer = new Timer(1000 / FPS, new ActionListener() {

			private int i = 2; //Skip some frames to give a head start

			private long startTime = System.currentTimeMillis();
			private int ANIMATION_FRAMES = 15;

			public void actionPerformed(ActionEvent arg0) {

				boolean skipFrame = (i < (ANIMATION_FRAMES - 1)) && System.currentTimeMillis() > startTime + (1000 / FPS) * i;
				boolean done = false;

				do {
					
					done = i >= ANIMATION_FRAMES;
					
					if (!done) {
						zoom(centerX, wheelRotation, skipFrame);
						i++;

					} else {
						stopAnimation();
					}
					
				} while (skipFrame && !done);
			}
		});

		mouseAnimationTimer.setRepeats(true);
		mouseAnimationTimer.setCoalesce(true);
		mouseAnimationTimer.start();
	}

	private void stopAnimation() {
		if (mouseAnimationTimer != null) {
			mouseAnimationTimer.stop();
			mouseAnimationTimer = null;
		}
	}

	protected void zoom(int lockedX, double wheelRotation, boolean disableDrawing) {

		// not all views are zoomed (e.g., the overview with cytoband)
		if (zoomable) {

			BpCoordDouble pointerBp = trackToBp(lockedX);
			double pointerRelative = trackToRelative(lockedX);

			double startBp = getBpRegionDouble().start.bp;
			double endBp = getBpRegionDouble().end.bp;

			double width = endBp - startBp;
			width *= Math.pow(ZOOM_FACTOR, wheelRotation);

			int minBpWidth = (int) (((float) parentPlot.chartPanel.getPreferredSize().getSize().width) / MIN_PIXELS_PER_NUCLEOTIDE);
			if (width < minBpWidth) {
				width = minBpWidth;
			}

			startBp = (double) (pointerBp.bp - width * pointerRelative);
			endBp = (double) (pointerBp.bp + width * (1 - pointerRelative));
			
			setBpRegion(new RegionDouble(startBp, getBpRegionDouble().start.chr, endBp, getBpRegionDouble().end.chr), disableDrawing);
			
			if (!disableDrawing) {
				parentPlot.redraw();
			}
		}
	}

	public Integer bpToTrack(BpCoord bp) {
		if (bpRegion.start.chr.equals(bp.chr)) {
			return (int) Math.round(((bp.bp - getBpRegionDouble().start.bp) * bpWidth()) + getX());

		} else {
			return null;
		}
	}

	/**
	 * Precisely convert bp coordinate to pixel position in this view. Rounding should be performed just before drawing.
	 * 
	 * @param bp
	 * @return
	 */
	public Float bpToTrackFloat(BpCoord bp) {
		if (bpRegion.start.chr.equals(bp.chr)) {
			return (float) ((bp.bp - getBpRegionDouble().start.bp) * bpWidth()) + getX();

		} else {
			return null;
		}
	}

	/**
	 * Calculates width of a single bp in pixels for this view. Number is a float, so the rounding should be performed just before drawing.
	 * 
	 * @return width of a single bp in pixels for this view.
	 */
	public Float bpWidth() {
		return getWidth() / (float) getBpRegionDouble().getLength();
	}

	public BpCoordDouble trackToBp(double d) {
		return new BpCoordDouble((double) (trackToRelative(d) * (getBpRegionDouble().end.bp - getBpRegionDouble().start.bp) + getBpRegionDouble().start.bp), bpRegion.start.chr);
	}

	public double trackToRelative(double track) {
		return (double) (track - getX()) / getWidth();
	}

	public int getX() {
		return viewArea.x;
	}

	public int getY() {
		return viewArea.y;
	}

	public void redraw() {
		//Dont accept redraw request from tracks if animation is running
		if (mouseAnimationTimer == null || !mouseAnimationTimer.isRunning()) {
			parentPlot.redraw();
		}
	}

	public void addRegionListener(RegionListener listener) {
		listeners.add(listener);
	}

	public void dispatchRegionChange() {
		for (RegionListener listener : listeners) {
			listener.regionChanged(getBpRegion());
		}
	}

	private Point2D scale(Point2D p) {
		return new Point((int) (p.getX() / parentPlot.chartPanel.getScaleX()), (int) (p.getY() / parentPlot.chartPanel.getScaleY()));
	}

	public String tooltipRequest(MouseEvent mouseEvent) {
		Point locationOnPanel = (Point) mouseEvent.getLocationOnScreen().clone();
		SwingUtilities.convertPointFromScreen(locationOnPanel, parentPlot.chartPanel);
		return tooltipRequest(locationOnPanel);
	}

	public String tooltipRequest(Point2D locationOnPanel) {
		return null; // tooltips disabled by default in views
	}

	public void clean() {
		
		if (queueManager != null) {
			queueManager.poisonAll();
		}
		scrollGroups.clear();
		//Queue manager holds references to track data through the data listener references preventing gc
		queueManager = null;
	}

	public void setViewLimiter(ViewLimiter viewLimiter) {

		this.viewLimiter = viewLimiter;
		viewLimiter.addLimitChangeListener(new RegionListener() {
			
			@Override
			public void regionChanged(Region bpRegion) {
				
				limitRegion();
				redraw();
			}

		});
	}
	
	public void limitRegion() {
		
		long maxBp = -1;
		boolean isEndLimited = false;
		
		if (viewLimiter != null && viewLimiter.getLimit() != null) {
			
			BpCoord limit = viewLimiter.getLimit();
			
			if (limit.chr.equals(bpRegion.start.chr)) {
				maxBp = getMaxBp(limit.bp);
			}
			
		} 
		
		if (maxBp == -1) {
			
			maxBp = (long)(double)this.bpRegion.end.bp;
		}
				
		RegionDouble limitedRegion = this.bpRegion.clone();
		
		if (limitedRegion.end.bp > maxBp) {
			limitedRegion.move(maxBp - limitedRegion.end.bp);
			isEndLimited = true;
		}
								
		//Enable scrolling to minus coordinates for 1/30 of width to 
		//make it easier to navigate to the beginning of chromosome  
		long minBp = getMinBp((long) Math.min(maxBp,this.bpRegion.end.bp));
		
		if (limitedRegion.start.bp < minBp ) {
			
			if (isEndLimited) {
				limitedRegion.start.bp = (double) minBp;				
			} else {
				limitedRegion.move(minBp-limitedRegion.start.bp);
			}
		}		
		
		this.bpRegion = limitedRegion;
	}
	
	public ViewLimiter getViewLimiter() {
		return this.viewLimiter; 
	}

	@Override
	public void setHeight(int height) {
		this.layoutHeight = height;
	}
	
	@Override
	public int getMinHeight() {
		return LayoutTool.getMinHeightSum(this);
	}
	
	@Override
	public boolean isVisible() {
		return true;
	}

	public void addScrollGroup(ScrollGroup group) {
		scrollGroups.add(group);
	}

	public Collection<? extends ScrollGroup> getScrollGroups() {
		return scrollGroups;
	}
	
	@Override
	public Collection<? extends LayoutComponent> getLayoutComponents() {
		return getScrollGroups();
	}
}
