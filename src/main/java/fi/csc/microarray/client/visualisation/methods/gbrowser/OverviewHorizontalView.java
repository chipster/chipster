package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.event.MouseEvent;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * Special version of horizontal view of tracks. It is used for the chromosome overview panel.
 *  
 * @author Petri Klemelä
 *
 */
public class OverviewHorizontalView extends HorizontalView {
	
	private List<RegionListener> overviewListeners = new LinkedList<RegionListener>();
	
	public OverviewHorizontalView(GenomePlot parent) {
		super(parent, false, false, false);
	}
	
	public void mouseClicked(MouseEvent e) {
		
		if (!viewArea.contains(e.getPoint())) {
			return;
		}
		
		Region newRegion = highlight;
		long delta = (long) (newRegion.getMid() - trackToBp(e.getX()).bp);
		newRegion.start.bp -= delta;
		newRegion.end.bp -= delta;
		
		dispatchOverviewRegionChange(newRegion);
		
	}

	public void mouseEntered(MouseEvent e) {
		// ignore
	}

	public void mouseExited(MouseEvent e) {
		// ignore
	}

	public void mousePressed(MouseEvent e) {

//		stopAnimation();
//		dragStartPoint = scale(e.getPoint());
//		dragStarted = false;
	}

	public void mouseReleased(MouseEvent e) {

//		if (dragStarted && dragEndPoint != null && dragLastStartPoint != null && Math.abs(dragEndPoint.getX() - dragLastStartPoint.getX()) > 10 && System.currentTimeMillis() - dragEventTime < DRAG_EXPIRATION_TIME_MS) {
//
//			stopAnimation();
//
//			timer = new Timer(1000 / FPS, new ActionListener() {
//
//				private int i = 0;
//				private int ANIMATION_FRAMES = 30;
//				private long startTime = System.currentTimeMillis();
//
//				public void actionPerformed(ActionEvent arg0) {
//
//					double endX = dragEndPoint.getX();
//					double startX = dragLastStartPoint.getX();
//
//					double newX = endX - (endX - startX) / (ANIMATION_FRAMES - i);
//
//					dragEndPoint = new Point2D.Double(newX, dragEndPoint.getY());
//
//					boolean skipFrame = (i < (ANIMATION_FRAMES - 1)) && System.currentTimeMillis() > startTime + (1000 / FPS) * i;
//
//					if (i < ANIMATION_FRAMES) {
//						handleDrag(dragLastStartPoint, dragEndPoint, skipFrame);
//						i++;
//					} else {
//						stopAnimation();
//					}
//				}
//			});
//			timer.setRepeats(true);
//			timer.start();
//		}
	}

	public void mouseDragged(MouseEvent e) {
		
//		if (movable) {
//			dragStarted = true;
//			dragEndPoint = scale(e.getPoint());
//			dragEventTime = System.currentTimeMillis();
//
//			handleDrag(dragStartPoint, dragEndPoint, false);
//
//		}
//		dragLastStartPoint = dragStartPoint;
//		dragStartPoint = scale(e.getPoint());
	}
	
	
	public void addOverviewRegionListener(RegionListener listener) {
		overviewListeners.add(listener);
	}

	public void dispatchOverviewRegionChange(Region selectedRegion) {
		for (RegionListener listener : overviewListeners) {
			listener.regionChanged(selectedRegion);
		}
	}
	
	@Override
	public boolean isRulerEnabled() {
		return false;
	}
}
