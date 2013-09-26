package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

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
	
	public OverviewHorizontalView(GBrowserPlot parent) {
		super(parent, false, false, false);
	}
	
	public void mouseClicked(MouseEvent e) {	
		
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
		// ignore
	}

	public void mouseReleased(MouseEvent e) {
		// ignore
	}

	public void mouseDragged(MouseEvent e) {
		// ignore		
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
	public boolean isCursorLineEnabled() {
		return false;
	}
}
