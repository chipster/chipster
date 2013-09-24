package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.Graphics;
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

import javax.swing.JPanel;
import javax.swing.Timer;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;

/**
 * Combines tracks to create a single synchronised view. All tracks within one view move horizontally hand-in-hand. 
 * View is responsible for handling user input that moves or zooms the view. 
 * Also the drawing of drawables is implemented in View, because different views 
 * (e.g. Horizontal or Circular) require different drawing implementations.
 * 
 */
public abstract class GBrowserView implements MouseListener, MouseMotionListener, MouseWheelListener { //, TooltipRequestProcessor {
	
	JPanel component = new JPanel() {			

		@Override
		public void paintChildren(Graphics g) {								
			
			super.paintChildren(g);
			
			// draw simple vertical cursor line
			if (isCursorLineEnabled()) {
				g.setColor(new Color(0f, 0f, 0f, 0.25f));
				g.drawLine((int)getWidth()/2, 0, (int)getWidth()/2, (int)getHeight());
			}
			
			// Show current position on top of chromosome cytoband
			if (highlight != null) {
				Rectangle rect = new Rectangle(g.getClip().getBounds());

				rect.x = bpToTrack(highlight.start);
				rect.width = Math.max(3, bpToTrack(highlight.end) - rect.x);
				rect.height = 20;

				g.setColor(new Color(0, 0, 0, 64));
				g.fillRect(rect.x, rect.y, rect.width, rect.height);
				g.setColor(new Color(0, 0, 0, 255));
				g.drawRect(rect.x, rect.y, rect.width, rect.height);
			}
		}
	};

	public RegionDouble bpRegion;
	public Region highlight;

	public Collection<ScrollGroup> scrollGroups = new LinkedList<ScrollGroup>();

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

	private Point2D dragEndPoint;
	private Point2D dragLastStartPoint;
	private long dragEventTime;
	private Region requestRegion;
	private Collection<Track> previousTracks;
	
	private Timer mouseAnimationTimer;

	private ViewLimiter viewLimiter;

	private static final long DRAG_EXPIRATION_TIME_MS = 50;

	public GBrowserView(GBrowserPlot parent, boolean movable, boolean zoomable, boolean selectable) {
		this.parentPlot = parent;
		this.movable = movable;
		this.zoomable = zoomable;	
				
		component.setLayout(new MigLayout("flowy, fillx, gap 0! 0!, insets 0"));
		component.setInheritsPopupMenu(true);
	}

	/**
	 * Add a track group containing one or several tracks.
	 * 
	 * @param group
	 */
	public void addTrackGroup(TrackGroup group) {
		ScrollGroup scrollGroup = new ScrollGroup();
		scrollGroup.addTrackGroup(group);
		addScrollGroup(scrollGroup);
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
				tracks.addAll(trackGroup.getTracks());
			}
		}
		return tracks;
	}

	public boolean isCursorLineEnabled() {
		return true;
	}

	public QueueManager getQueueManager() {	
		if (queueManager == null) {
			queueManager = new QueueManager(this);
		}
		return queueManager;
	}

	/**
	 * Fire data requests for all tracks in this view. Request is made only if the previous request
	 * doesn't cover current area. Requests are enlarged to avoid an immediate need for new request when
	 * the view is moved. However, requests can not be enlarged too much, because the result data has to
	 * fit in the memory.
	 * 
	 * Data is requested from the region 
	 * which is twice the size of the visible region, i.e. it extends to the region outside 
	 * visible region on both sides for the length of half of the view size. A new request 
	 * is made when the previous request covers only a quarter of the view length outside 
	 * the visible region. This gives us little bit time to receive and process the data before the user
	 * has moved outisde the data area of the previous request.
	 * 
	 * Only fire one request for a single file.
	 */
	public void fireDataRequests() {

		Map<DataThread, Set<DataType>> datas = new HashMap<DataThread, Set<DataType>>();

		// Add all requested columns for each requested file
		for (Track t : getTracks()) { 
					
			t.clearDataTypes();
			t.defineDataTypes();
			
			Map<DataThread, Set<DataType>> trackDatas = t.getDataTypeMap();
			
			boolean cancelType = false;
			
			// Don't do anything for tracks without data
			if (trackDatas == null) {
				continue;
			} 
				
			for (Set<DataType> types : trackDatas.values()) {
				if (types.contains(DataType.CANCEL)) {
					cancelType = true;
				}
			}

			// Don't do anything for hidden tracks, unless they wan't to send cancel request
			if (!cancelType && !t.isSuitableViewLength()) {
				continue;
			}

			for (DataThread file : trackDatas.keySet()) {
				
				if (file != null) {
					// Add columns for this requested file
					Set<DataType> columns = datas.get(file);
					columns = columns != null ? columns : new HashSet<DataType>();
					columns.addAll(trackDatas.get(file));
					datas.put(file, columns);
				}
			}
		}
		
		Region view = getBpRegion();
		
		long minBuffer = view.getLength() / 4;
		long maxBuffer = view.getLength() / 2;
		
		long minStart = view.start.bp - maxBuffer;
		long maxStart = view.start.bp - minBuffer;	
		long minEnd = view.end.bp + minBuffer;
		long maxEnd = view.end.bp + maxBuffer;
		
		boolean unchangedTracks = unchangedTracks();
		
		if (unchangedTracks && requestRegion != null && view.start.chr.equals(requestRegion.start.chr) &&
				requestRegion.start.bp < maxStart && requestRegion.end.bp > minEnd &&
				getBpRegion().getLength() > requestRegion.getLength() / 4) {
			
			//Previous request is enough
			
		} else {
			
			//Get new data
			
			Region newRequest = new Region(minStart, maxEnd, view.start.chr);
			Region previousRequest = requestRegion;
			requestRegion = new Region(newRequest.start, newRequest.end);
						
			//if we have requested part of this data already, don't request that part again
			if (unchangedTracks && previousRequest != null && newRequest.start.chr.equals(previousRequest.start.chr) && newRequest.intersects(previousRequest)) {
				
				Region overlap = newRequest.intersect(previousRequest);

				if (overlap.start.equals(newRequest.start)) {
					// Overlaps from left
					newRequest.start = overlap.end;

				} else if (overlap.end.equals(newRequest.end)) {
					// Overlaps from right
					newRequest.end = overlap.start;

				} else {
					// Overlap inside request, do nothing, because would need splitting
				}
			} 
						
			// Fire data requests
			for (DataThread file : datas.keySet()) {
				DataStatus status = new DataStatus();
				getQueueManager().addDataRequest(file, new DataRequest(newRequest, datas.get(file), status), requestRegion);
			}
		}
	}
	
	private boolean unchangedTracks() {
		
		Collection<Track> currentTracks = getVisibleTracks();
		
		boolean unchanged = false;
		
		if (previousTracks !=  null && 
				previousTracks.containsAll(currentTracks) &&
					currentTracks.containsAll(previousTracks)) {
			
			unchanged = true;			
		} else {
			
			unchanged = false;
		}	
		
		previousTracks = getVisibleTracks();		
		return unchanged;
	}

	private LinkedList<Track> getVisibleTracks() {
		
		LinkedList<Track> list = new LinkedList<Track>();
		
		for (Track track : getTracks()) {
			if (track.isSuitableViewLength()) {
				list.add(track);
			}
		}
		
		return list;
	}

	public long getMinBp(long length) {
		return (long) (-length / 30);
	}
	
	public long getMaxBp(long length) {
		return (long) (length * (1 + 1.0 / 30));
	}

	public void setBpRegion(RegionDouble region) {				

		this.bpRegion = limitRegion(region);		
		
		fireDataRequests();
		
		dispatchRegionChange();
	}

	protected void updateLayout() {
		
		component.removeAll();
				
		for (ScrollGroup group : scrollGroups) {
			group.updateLayout();
			
	        LayoutMode mode = group.getLayoutMode();
	        if (LayoutMode.FIXED == mode) {	        	
	        	component.add(group.getComponent(), "growx, ");
	        } else {
	        	component.add(group.getComponent(), "grow");
	        }		
		}
		
//		printLayout();
		
		//ScrollPane sizes may have changed
		component.revalidate();
	}

	
	/**
	 * Print layout settings for debugging
	 */
//	private void printLayout() {
//		for (ScrollGroup group : scrollGroups) {
//			if ("Samples".equals(group.getScrollGroupName()) || "Annotations".equals(group.getScrollGroupName())) {
//				System.out.println(
//						group.getScrollGroupName() + 
//						"\tgetSize: " + group.getComponent().getHeight() +  
//						"\tgetPreferredSize: " + group.getComponent().getPreferredSize().getHeight() + 
//						"\tgetLayoutMode: " + group.getLayoutMode());
//
//				for (TrackGroup trackGroup : group.getTrackGroups()) {
//					System.out.println(
//							"\t" + trackGroup.getComponent().getName() + 
//							"\tgetSize: " + trackGroup.getComponent().getSize().getHeight() + 
//							"\tgetPreferredSize: " + trackGroup.getComponent().getPreferredSize().getHeight() + 
//							"\tgetLayoutMode: " + trackGroup.getLayoutMode());
//					
//					for (Track track : trackGroup.getTracks()) {
//						System.out.println(
//								"\t\t" + track.getComponent().getName() + 
//								"\tgetSize: " + track.getComponent().getSize().getHeight() + 
//								"\tgetPreferredSize: " + track.getComponent().getPreferredSize().getHeight() + 
//								"\tgetLayoutMode: " + track.getLayoutMode());					
//					}
//				}
//			}
//		}
//	}

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

		component.requestFocusInWindow();

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
		
		if (movable && dragStartPoint != null) {

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

			int minBpWidth = (int) (getWidth() / MIN_PIXELS_PER_NUCLEOTIDE);
			if (width < minBpWidth) {
				width = minBpWidth;
			}

			startBp = (double) (pointerBp.bp - width * pointerRelative);
			endBp = (double) (pointerBp.bp + width * (1 - pointerRelative));
			
			setBpRegion(new RegionDouble(startBp, getBpRegionDouble().start.chr, endBp, getBpRegionDouble().end.chr));
			
			if (!disableDrawing) {
				parentPlot.redraw();
			}
		}
	}

	public Integer bpToTrack(BpCoord bp) {
		if (bpRegion.start.chr.equals(bp.chr)) {
			return (int) Math.round(((bp.bp - getBpRegionDouble().start.bp) * bpWidth()));

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
			return (float) ((bp.bp - getBpRegionDouble().start.bp) * bpWidth());

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
		return (double) (track) / getWidth();
	}

	public int getWidth() {
		return component.getWidth();
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
		/*Scaling is used only when the visualization is printed or saved as an image.
		 * In that case, the Java Graphics library does the scaling automatically. In interactive mode
		 * the component is rendered always without scaling, so there is not need for scaling the
		 * coordinates of the user input at the moment.   
		 */
		return p;
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
				
				setBpRegion(getBpRegionDouble());
				redraw();
			}

		});
	}
	
	public RegionDouble limitRegion(RegionDouble region) {
		
		long maxBp = -1;
		boolean isEndLimited = false;
		
		if (viewLimiter != null && viewLimiter.getLimit() != null) {
			
			BpCoord limit = viewLimiter.getLimit();
			
			if (limit.chr.equals(region.start.chr)) {
				maxBp = getMaxBp(limit.bp);
			}
			
		} 
		
		if (maxBp == -1) {
			
			maxBp = (long)(double)region.end.bp;
		}
				
		RegionDouble limitedRegion = region.clone();
		
		if (limitedRegion.end.bp > maxBp) {
			limitedRegion.move(maxBp - limitedRegion.end.bp);
			isEndLimited = true;
		}
								
		//Enable scrolling to minus coordinates for 1/30 of width to 
		//make it easier to navigate to the beginning of chromosome  
		long minBp = getMinBp((long) Math.min(maxBp,region.end.bp));
		
		if (limitedRegion.start.bp < minBp ) {
			
			if (isEndLimited) {
				limitedRegion.start.bp = (double) minBp;				
			} else {
				limitedRegion.move(minBp-limitedRegion.start.bp);
			}
		}		
		
		return limitedRegion;
	}
	
	public ViewLimiter getViewLimiter() {
		return this.viewLimiter; 
	}

	public void addScrollGroup(ScrollGroup group) {
		scrollGroups.add(group);
		group.setView(this);
		this.updateLayout();
	}

	public Collection<? extends ScrollGroup> getScrollGroups() {
		return scrollGroups;
	}

	/**
	 * Use this method to check if a data item should be kept or removed. Data requests aren't equal to visible 
	 * region. The data items returned by the last request won't be requested again and therefore those must not be
	 * removed even if those aren't currently visible.   
	 * 
	 * @param region of the data item
	 * @return true if data item with the given region is still relevant or false if can be removed
	 */
	public boolean requestIntersects(Region region) {
		if (requestRegion != null) {			
			return requestRegion.intersects(region);
		} 
		
		//Keep everything
		return true;
	}

	/**
	 * Similar to @link GBrowserView#requestIntersects(Region), but works with single position instead of region.
	 * 
	 * @param position
	 * @return
	 */
	public boolean requestContains(BpCoord position) {
		if (requestRegion != null) {			
			return requestRegion.contains(position);
		} 
		
		//Keep everything
		return true;		
	}
	
	/**
	 * For checking if the data is still needed, use @link GBrowserView#requestIntersects(Region).
	 * 
	 * @param region
	 * @return
	 */
	public boolean viewIntersects(Region region) {
					
		return getBpRegion().intersects(region);		
	}


	public Region getRequestRegion() {
		return requestRegion;
	}

	public void reloadData() {
		getQueueManager().clearDataResultListeners();
		parentPlot.getBrowser().initializeDataResultListeners();
		reloadDataLater();
		fireDataRequests();
		redraw();
	}

	public void reloadDataLater() {
		requestRegion = null;
	}
	
	public LayoutMode getLayoutMode() {
		return LayoutTool.inferViewLayoutMode(scrollGroups);
	}

	public JPanel getComponent() {
		return component;
	}

	public GBrowserPlot getPlot() {
		return parentPlot;
	}

	public LinkedList<TrackGroup> getTrackGroups() {
		LinkedList<TrackGroup> groups = new LinkedList<>();
		for (ScrollGroup scrollGroup : getScrollGroups()) {
			groups.addAll(scrollGroup.getTrackGroups());
		}
		return groups;
	}
}
