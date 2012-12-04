package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Insets;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
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

	public class ScrollGroupBar extends JScrollBar implements AdjustmentListener {
		
		public ScrollGroupBar() {
			addAdjustmentListener(this);
		}
		
		private int referenceY = 0;
			
		public void setHeight(int height) {
			int scrollBarWidth = ((Integer) UIManager.get("ScrollBar.width")).intValue();
			setSize(scrollBarWidth, height);
		}
		
		/**
		 * Set maximum, but keep the same scroll position in comparison to referenceY
		 * 
		 * @param maximum
		 * @param extent
		 * @param referenceY
		 */
		public void set(int maximum, int extent, int referenceY) {
			int offset = getValue() - referenceY + getModel().getExtent() / 2;
			setMaximum(maximum);
			getModel().setExtent(extent);
			setValue(referenceY - extent / 2 + offset);
			this.referenceY = referenceY;
		}
		
		public void setDefaultValue() {
			setValue(referenceY - getModel().getExtent() / 2);
		}
		
		@Override
		public void adjustmentValueChanged(AdjustmentEvent e) {
			GBrowserChartPanel.this.genomePlot.redraw();
		}
	}

	public GBrowserChartPanel() {
		super(null);		
	}

	private GBrowserPlot genomePlot;
	//FIXME clean when visualization is closed
	private Map<ScrollGroup, ScrollGroupBar> scrollBarsMap = new HashMap<ScrollGroup, ScrollGroupBar>();
	//Preserves ScrollGroup order
	private List<ScrollGroupBar> scrollBarsList = new LinkedList<ScrollGroupBar>();

	public void paintComponent(Graphics g) {

		//in genomeBrowser getSize() doesn't work below width of 680
		//but parent seems to know better

		int height =  getParent().getSize().height;

		/* If the genomeBrowser height is bigger than the scrollPane height, we set that bigger value
		 * to both chartpanel size and drawHeight to make the genomeBrowser draw vertically everything 
		 * without scaling
		 */
//		if (genomePlot != null && genomePlot.isLegacyFullHeight()) {
//
//			height = genomePlot.getHeight();
//		}

//		if (getParent() instanceof JViewport) {
//			JViewport viewport = (JViewport)getParent();
//			genomePlot.setLegacyFullHeightClip(viewport.getViewRect());
//		}

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

			//Calculate new values, but don't set them yet, because we still need the old values
			int barHeight = Math.min(group.getHeight(), maxY - barHeightSum);
			barHeightSum += barHeight;
			int extent = barHeight;
			
			int maximum;
			if (group.isScrollEnabled()) {
				maximum = group.getFullHeight();
			} else {
				maximum = group.getHeight();
			}

			//Search for the scroll bar of this ScrollGroup, create a new one if it doesn't exist yet 
			ScrollGroupBar bar;
			if (!scrollBarsMap.containsKey(group)) {
				bar = new ScrollGroupBar();
				
				scrollBarsMap.put(group, bar);

			} else {					
				bar = scrollBarsMap.get(group); 
			}

			boolean visible = group.isVisible() && group.isScrollEnabled() && maximum > barHeight;
			boolean becomesVisible = bar.isVisible() == false && visible == true;
			
			int referenceY = group.getScrollReferenceY();

			//FIXME replace with referenceY
			if (ScrollPosition.START == group.getDefaultScrollPosition()) {
				if (becomesVisible) {
					bar.setValue(0);
				}
				bar.set(maximum, extent, referenceY);
			}

			if (ScrollPosition.MID == group.getDefaultScrollPosition()) {
				if (becomesVisible) {
					bar.setDefaultValue();
				}
				bar.set(maximum, extent, referenceY);
			}

			//Set the new values
			bar.setHeight(barHeight);
			bar.setVisible(visible);

			scrollBarsList.add(bar);			
		}
	}

	public int getScrollValue(ScrollGroup scrollGroup) {
		if (scrollBarsMap.containsKey(scrollGroup)) {
			return scrollBarsMap.get(scrollGroup).getValue();
		} else {
			return 0;
		}
	}
}
