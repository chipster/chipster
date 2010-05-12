package fi.csc.microarray.client.visualisation.methods.gbrowser;

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
import java.awt.image.BufferedImage;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import javax.swing.Timer;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.QueueManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.RulerTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;

public abstract class View implements MouseListener, MouseMotionListener, MouseWheelListener {

	protected BpCoordRegionDouble bpRegion;
	public BpCoordRegion highlight;

	public Collection<Track> tracks = new LinkedList<Track>();
	protected Rectangle viewArea = new Rectangle(0, 0, 500, 500);
	private QueueManager queueManager = new QueueManager();
	private Point2D dragStartPoint;
	private boolean dragStarted;

	GenomePlot parentPlot;

	private static final int FPS = 30;

	private boolean movable;
	protected boolean zoomable;

	protected final float ZOOM_FACTOR = 1.06f;

	private List<RegionListener> listeners = new LinkedList<RegionListener>();
	public int margin = 0;
	protected Float trackHeight;
	private Point2D dragEndPoint;
	private Point2D dragLastStartPoint;
	private Iterator<Track> trackIter;
	private Iterator<Drawable> drawableIter;
	private boolean continueDrawingLater;
	private Track track;
	private float y;
	private int x;
	private BufferedImage drawBuffer;
	private long dragEventTime;
	private static final long DRAG_EXPIRATION_TIME_MS = 50;

	public View(GenomePlot parent, boolean movable, boolean zoomable, boolean selectable) {
		parentPlot = parent;
		this.movable = movable;
		this.zoomable = zoomable;
	}

	protected abstract void drawDrawable(Graphics2D g, int x, int y, Drawable drawable);


	public void addTrack(Track t) {
		tracks.add(t);
	}

	protected void drawView(Graphics2D g, boolean isAnimation) {

		if (bpRegion == null) {
			setBpRegion(new BpCoordRegionDouble(0d, 1024 * 1024 * 250d, new Chromosome("1")), false);
		}

		Rectangle viewClip = g.getClipBounds();
		viewArea = viewClip;

		if (drawBuffer == null) {
			drawBuffer = new BufferedImage((int) viewArea.getWidth(), (int) viewArea.getHeight(), BufferedImage.TYPE_INT_RGB);

			Graphics2D bufG2 = (Graphics2D) drawBuffer.getGraphics();
			bufG2.setPaint(Color.white);
			bufG2.fillRect(0, 0, drawBuffer.getWidth(), drawBuffer.getHeight());
		}

		Graphics2D bufG2 = (Graphics2D) drawBuffer.getGraphics();

		if (trackIter == null) {

			// y = viewClip.y;
			// x = viewClip.x;

			y = 0;
			x = 0;

			trackIter = tracks.iterator();
			drawableIter = null;
		}

//		long startTime = System.currentTimeMillis();
		continueDrawingLater = false;

		while (trackIter.hasNext() || (drawableIter != null && drawableIter.hasNext())) {

			if (drawableIter == null || !drawableIter.hasNext()) {
				track = trackIter.next();
			}

			if (track.getMaxHeight() > 0) {

				if (drawableIter == null) {
					Collection<Drawable> drawables = track.getDrawables();
					drawableIter = drawables.iterator();
				}

				while (drawableIter.hasNext()) {

					Drawable drawable = drawableIter.next();
					
					if(drawable == null) {
						continue;
					}

					int maybeReversedY = (int) y;

					if (track.getStrand() == Strand.REVERSED) {
						drawable.upsideDown();
						maybeReversedY += Math.min(getTrackHeight(), track.getMaxHeight());
					}

					drawDrawable(bufG2, x, maybeReversedY, drawable);

//					if (System.currentTimeMillis() - startTime >= 1000 / FPS) {
//						continueDrawingLater = true;
//						this.redraw();
//						break;
//					}
				}

				if (continueDrawingLater) {
					break;
					
				} else {
					drawableIter = null;
				}
				
			} else {
				drawableIter = null;
			}

			if (track.getMaxHeight() == Integer.MAX_VALUE) {
				y += getTrackHeight();
			} else {
				y += track.getMaxHeight();
			}
		}

		g.drawImage(drawBuffer, (int) viewArea.getX(), (int) viewArea.getY(), (int) viewArea.getX() + drawBuffer.getWidth(), (int) viewArea.getY() + drawBuffer.getHeight(), 0, 0, drawBuffer.getWidth(), drawBuffer.getHeight(), null);

		if (!continueDrawingLater) {
			bufG2.setPaint(Color.white);
			bufG2.fillRect(0, 0, drawBuffer.getWidth(), drawBuffer.getHeight());
			trackIter = null;
			drawableIter = null;
		}

	}

	public float getTrackHeight() {

		if (trackHeight == null) {
			trackHeight = (getHeight() - getMaxTrackHeightTotal()) / (float) getStretchableTrackCount();
		}
		
		return trackHeight;
	}

	protected int getMaxTrackHeightTotal() {
		int maxHeightTotal = 0;

		for (Track track : tracks) {
			if (track.getMaxHeight() != Integer.MAX_VALUE) {
				maxHeightTotal += track.getMaxHeight();
			}
		}
		return maxHeightTotal;
	}

	protected int getStretchableTrackCount() {
		int stretchableCount = 0;

		for (Track track : tracks) {
			if (track.getMaxHeight() == Integer.MAX_VALUE) {
				stretchableCount++;
			}
		}
		return stretchableCount;
	}

	public int getWidth() {
		return this.viewArea.width;
	}

	public int getHeight() {
		return this.viewArea.height;
	}

	public QueueManager getQueueManager() {
		return queueManager;
	}

	public void setBpRegion(BpCoordRegionDouble region, boolean disableDrawing) {
		bpRegion = region;

		// Bp-region change may change visibility of tracks, calculate sizes again
		trackHeight = null;

		if (!disableDrawing) {
			for (Track t : tracks) {

				t.updateData();
			}
			dispatchRegionChange();
		}
	}

	public BpCoordRegionDouble getBpRegionDouble() {
		return bpRegion;
	}

	public BpCoordRegion getBpRegion() {
		return new BpCoordRegion((long) (double) bpRegion.start.bp, bpRegion.start.chr, (long) (double) bpRegion.end.bp, bpRegion.end.chr);
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

		stopAnimation();
		dragStartPoint = scale(e.getPoint());
		dragStarted = false;
	}

	public void mouseReleased(MouseEvent e) {

		if (dragStarted && dragEndPoint != null && dragLastStartPoint != null && Math.abs(dragEndPoint.getX() - dragLastStartPoint.getX()) > 10 && System.currentTimeMillis() - dragEventTime < DRAG_EXPIRATION_TIME_MS) {

			stopAnimation();

			timer = new Timer(1000 / FPS, new ActionListener() {

				private int i = 0;
				private int ANIMATION_FRAMES = 30;
				private long startTime = System.currentTimeMillis();

				public void actionPerformed(ActionEvent arg0) {

					double endX = dragEndPoint.getX();
					double startX = dragLastStartPoint.getX();

					double newX = endX - (endX - startX) / (ANIMATION_FRAMES - i);

					dragEndPoint = new Point2D.Double(newX, dragEndPoint.getY());

					boolean skipFrame = (i < (ANIMATION_FRAMES - 1)) && System.currentTimeMillis() > startTime + (1000 / FPS) * i;

					if (i < ANIMATION_FRAMES) {
						handleDrag(dragLastStartPoint, dragEndPoint, skipFrame);
						i++;
					} else {
						stopAnimation();
					}
				}
			});
			timer.setRepeats(true);
			timer.start();
		}
	}

	public void mouseDragged(MouseEvent e) {
		
		if (movable) {
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

	private Timer timer;

	public void mouseWheelMoved(final MouseWheelEvent e) {

		stopAnimation();

		timer = new Timer(1000 / FPS, new ActionListener() {

			private int i = 0;

			// 100 ms to give time for slower machines to view couple animation frames also
			private long startTime = System.currentTimeMillis() + 100;
			private int ANIMATION_FRAMES = 15;

			public void actionPerformed(ActionEvent arg0) {

				boolean skipFrame = (i < (ANIMATION_FRAMES - 1)) && System.currentTimeMillis() > startTime + (1000 / FPS) * i;

				if (i < ANIMATION_FRAMES) {
					zoom((int) scale(e.getPoint()).getX(), e.getWheelRotation(), skipFrame);
					i++;
					
				} else {
					stopAnimation();
				}
			}
		});
		
		timer.setRepeats(true);
		timer.setCoalesce(false);
		timer.start();
	}

	private void stopAnimation() {
		if (timer != null) {
			timer.stop();
			timer = null;
		}
	}

	private void zoom(int lockedX, int wheelRotation, boolean disableDrawing) {

		if (zoomable) {
			
			if (wheelRotation > 0) {
				lockedX = (int) getWidth() - lockedX + getX() * 2;
			}

			BpCoordDouble pointerBp = trackToBp(lockedX);
			double pointerRelative = trackToRelative(lockedX);

			double startBp = getBpRegionDouble().start.bp;
			double endBp = getBpRegionDouble().end.bp;

			double width = endBp - startBp;
			width *= Math.pow(ZOOM_FACTOR, wheelRotation);

			startBp = (double) (pointerBp.bp - width * pointerRelative);
			endBp = (double) (pointerBp.bp + width * (1 - pointerRelative));

			if (startBp < 0) {
				endBp += -startBp;
				startBp = 0;
			}

			setBpRegion(new BpCoordRegionDouble(startBp, getBpRegionDouble().start.chr, endBp, getBpRegionDouble().end.chr), disableDrawing);
		}
	}

	public Integer bpToTrack(BpCoord bp) {

		if (bpRegion.start.chr.equals(bp.chr)) {
			return (int) (((double) bp.bp - getBpRegionDouble().start.bp) / (getBpRegionDouble().end.bp - getBpRegionDouble().start.bp) * getWidth()) + getX();
			
		} else {
			return null;
		}
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
		parentPlot.redraw();
	}

	public List<Long> getRulerInfo() {
		for (Track t : tracks) {
			if (t instanceof RulerTrack) {
				RulerTrack ruler = (RulerTrack) t;
				return ruler.getRulerInfo();
			}
		}
		return null;
	}

	public void addRegionListener(RegionListener listener) {
		listeners.add(listener);
	}

	public void dispatchRegionChange() {
		for (RegionListener listener : listeners) {
			listener.RegionChanged(getBpRegion());
		}
	}

	private Point2D scale(Point2D p) {
		return new Point((int) (p.getX() / parentPlot.chartPanel.getScaleX()), (int) (p.getY() / parentPlot.chartPanel.getScaleY()));
	}
}
