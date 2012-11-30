package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.geom.Rectangle2D;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.swing.JScrollBar;
import javax.swing.JViewport;
import javax.swing.UIManager;

import org.jfree.chart.ChartPanel;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.ScrollGroup.ScrollPosition;

/**
 * ChartPanel without internal scaling to be used inside JScrollPane.
 * 
 * @author klemela
 */
public class GBrowserChartPanel extends ChartPanel {

	public GBrowserChartPanel() {
		super(null);		
	}

	private GBrowserPlot genomePlot;
	//FIXME memory leak or not?
	private Map<ScrollGroup, JScrollBar> scrollBarsMap = new HashMap<ScrollGroup, JScrollBar>();
	//Preserves ScrollGroup order
	private List<JScrollBar> scrollBarsList = new LinkedList<JScrollBar>();

	public void paintComponent(Graphics g) {

		//in genomeBrowser getSize() doesn't work below width of 680
		//but parent seems to know better

		int height =  getParent().getSize().height;

		/* If the genomeBrowser height is bigger than the scrollPane height, we set that bigger value
		 * to both chartpanel size and drawHeight to make the genomeBrowser draw vertically everything 
		 * without scaling
		 */
		if (genomePlot != null && genomePlot.isFullHeight()) {

			height = genomePlot.getHeight();
		}

		if (getParent() instanceof JViewport) {
			JViewport viewport = (JViewport)getParent();
			genomePlot.setFullHeightClip(viewport.getViewRect());
		}

		Dimension size = new Dimension(getParent().getSize().width, height);

		//Required for JScrollPane to understand that content requires scrolling
		this.setPreferredSize(size);

		//Required for JFreeChart to avoid scaling after the view is drawn
		this.setSize(size);

		//Required for JFreeChart to do the view drawing in correct resolution
		Insets insets = getInsets();
		Rectangle2D available = new Rectangle2D.Double(insets.left, insets.top,
				size.getWidth() - insets.left - insets.right,
				size.getHeight() - insets.top - insets.bottom);

		this.setMinimumDrawWidth((int)available.getWidth() - 1);
		this.setMinimumDrawHeight((int)available.getHeight() - 1);
		this.setMaximumDrawWidth((int)available.getWidth() + 1);
		this.setMaximumDrawHeight((int)available.getHeight() + 1);

		udpateScrollBars();

		super.paintComponent(g);
	}

	private void udpateScrollBars() {

		this.removeAll();

		int maxY = getHeight() - getInsets().bottom;
		int scrollBarWidth = ((Integer) UIManager.get("ScrollBar.width")).intValue();
		int scrollBarX = getWidth() - getInsets().right - scrollBarWidth;

		int y = getInsets().top;

		for (JScrollBar bar : scrollBarsList) {

			if (bar != null) {
				int height = Math.min(maxY, bar.getHeight());

				bar.setBounds(
						scrollBarX,
						y,
						scrollBarWidth,
						height);

				y += height;

				this.add(bar);
			}
		}
	}

	public void setGenomePlot(GBrowserPlot plot) {
		genomePlot = plot;		
	}	

	public void setScrollGroupBoundaries(Collection<ScrollGroup> scrollGroups, int maxY) {

		//Keep only scroll bars that still have a ScrollGroup
		Iterator<ScrollGroup> iter = scrollBarsMap.keySet().iterator();
		while (iter.hasNext()) {
			ScrollGroup mapGroup = iter.next();
			if(!scrollGroups.contains(mapGroup)) {
				iter.remove();
			}
		}
		scrollBarsList.clear();

		//Calculate height sum to prevent scrollbars from going outside the window
		int barHeightSum = 0;

		for ( ScrollGroup group : scrollGroups) {			

			//Calculate new values, but don't set them yet, because the old values are needed
			int barHeight = Math.min(group.getHeight(), maxY - barHeightSum);
			barHeightSum += barHeight;
			int extent = barHeight;
			
			int maximum;
			if (group.isScrollEnabled()) {
				maximum = group.getCanvasHeight();
			} else {
				maximum = group.getHeight();
			}

			//Search for the scroll bar of this ScrollGroup, create a new one if it doesn't exist yet 
			JScrollBar bar;
			if (!scrollBarsMap.containsKey(group)) {
				bar = new JScrollBar(JScrollBar.VERTICAL);
				scrollBarsMap.put(group, bar);

			} else {					
				bar = scrollBarsMap.get(group); 
			}

			boolean visible = group.isVisible() && group.isScrollEnabled() && bar.getMaximum() > bar.getHeight();
			boolean becomesVisible = bar.isVisible() == false && visible == true;
			
			int referenceY = group.getScrollReferenceY();

			//FIXME replace with reference
			if (ScrollPosition.START == group.getDefaultScrollPosition()) {
				if (becomesVisible) {
					bar.setValue(0);
				}
			}

			if (ScrollPosition.MID == group.getDefaultScrollPosition()) {
				if (becomesVisible) {
					bar.setValue(getMidValue(maximum, barHeight));
				} else {
					//FIXME last reference value needed to calculate offset
					//Try to keep the old scroll position despite the scroll area has changed
//					int offset = bar.getValue() - getMidValue(bar.getMaximum(), bar.getModel().getExtent());
//					bar.setMaximum(maximum);
//					bar.setValue(getMidValue(maximum, extent) + offset);

					bar.setMaximum(maximum);
					//FIXME not quite right
					//bar.setValue(referenceY - extent / 2);
					
				}
			}

			//Set the new values
			int scrollBarWidth = ((Integer) UIManager.get("ScrollBar.width")).intValue();
			bar.setSize(scrollBarWidth, barHeight);
			bar.getModel().setExtent(extent);
			bar.setVisible(visible);

			scrollBarsList.add(bar);			
		}
	}

	private int getMidValue(int maximum, int extent) {
		return (maximum - extent) / 2;
	}
	

	public int getScrollValue(ScrollGroup scrollGroup) {
		if (scrollBarsMap.containsKey(scrollGroup)) {
			return scrollBarsMap.get(scrollGroup).getValue();
		} else {
			return 0;
		}
	}
}
