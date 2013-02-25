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
import javax.swing.UIManager;

import org.jfree.chart.ChartPanel;

/**
 * ChartPanel without internal scaling to be used inside JScrollPane.
 * 
 * @author klemela
 */
public class GBrowserChartPanel extends ChartPanel {

	/**
	 * GUI component for ScrollGroup. The actual scrolling happens in ScrollGroup.draw(...) method.
	 * The ChartPanel is lowest Swing component in the genome browser and therefore the scroll bars
	 * have to be created here separately from the ScrollGroup.
	 * 
	 * @author klemela
	 */
	public class ScrollGroupBar extends JScrollBar implements AdjustmentListener {
		
		public ScrollGroupBar() {
			addAdjustmentListener(this);
			this.setUnitIncrement(10);
		}
		
		/**
		 * Reference level of the scroll content. The reference level is kept steady between frames. In practice
		 * this is used to keep the Ruler track in place, although the amount of Transcripts varies 
		 * around it.
		 */
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
			
			int referenceDelta = referenceY - this.referenceY;
			this.referenceY = referenceY;

			setMaximum(maximum);
			getModel().setExtent(extent);
			setValue(this.getValue() + referenceDelta);						
		}
		
		/**
		 * Set scrolling position according to referenceY and forget any user induced 
		 * value changes.
		 */
		public void setDefaultValue() {
			setValue(referenceY - getModel().getExtent() / 2);
		}
		
		/* 
		 * Redraw plot to show the new scrolling position
		 */
		@Override
		public void adjustmentValueChanged(AdjustmentEvent e) {
			GBrowserChartPanel.this.plot.redraw();
		}
	}

	public GBrowserChartPanel() {
		super(null);		
	}

	private GBrowserPlot plot;
	//FIXME clean when visualization is closed
	private Map<ScrollGroup, ScrollGroupBar> scrollBarsMap = new HashMap<ScrollGroup, ScrollGroupBar>();
	//Preserves ScrollGroup order
	private List<ScrollGroupBar> scrollBarsList = new LinkedList<ScrollGroupBar>();
	private int maxY;

	public void paintComponent(Graphics g) {

		//in genomeBrowser getSize() doesn't work below width of 680
		//but parent seems to know better
		int height =  getParent().getSize().height;

		Dimension size = new Dimension(getParent().getSize().width, height);

		//Required for JScrollPane (not used in the current Chipster implementation) to understand that content requires scrolling
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

	/**
	 * Put scroll bars to this Swing container. Scroll bars are configured 
	 * in method setScrollGroupBoundaries(). 
	 */
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

	public void setPlot(GBrowserPlot plot) {
		this.plot = plot;		
	}	
	
	/**
	 * Use ScrollGroup view port and content sizes to calculate scroll bar values.
	 * 
	 * @param scrollGroups
	 * @param maxY
	 */
	public void setScrollGroupOrder(Collection<ScrollGroup> scrollGroups, int maxY) {
		
		this.maxY = maxY;

		//Keep only scroll bars that still have a ScrollGroup
		Iterator<ScrollGroup> iter = scrollBarsMap.keySet().iterator();
		while (iter.hasNext()) {
			ScrollGroup mapGroup = iter.next();
			if(!scrollGroups.contains(mapGroup)) {
				iter.remove();
			}
		}
		scrollBarsList.clear();

		for ( ScrollGroup group : scrollGroups) {			

			//Search for the scroll bar of this ScrollGroup, create a new one if it doesn't exist yet 
			ScrollGroupBar bar;
			if (!scrollBarsMap.containsKey(group)) {
				bar = new ScrollGroupBar();
				
				scrollBarsMap.put(group, bar);
			} else {					
				bar = scrollBarsMap.get(group); 
			}
			scrollBarsList.add(bar);			
		}
	}

	/**
	 * Use ScrollGroup view port and content sizes to calculate scroll bar values.
	 * 
	 * @param scrollGroups
	 * @param maxY
	 */
	public void setScrollGroupBoundaries(ScrollGroup group) {

		int barHeight = group.getHeight();

		int extent = barHeight;

		int maximum;
		if (group.isScrollEnabled()) {
			maximum = group.getFullHeight();
		} else {
			maximum = group.getHeight();
		}

		ScrollGroupBar bar = scrollBarsMap.get(group);
		if (bar != null) {

			boolean visible = group.isVisible() && group.isScrollEnabled() && maximum > barHeight;			
			int referenceY = group.getScrollReferenceY();

			//Set the new values
			bar.set(maximum, extent, referenceY);
			bar.setHeight(barHeight);
			bar.setVisible(visible);

			this.validate();
		}
	}

	/**
	 * This method is used in ScrollGroup.draw() to get the scrolling position of each ScrollGroup
	 * to show the right part of the content.
	 * 
	 * @param scrollGroup
	 * @return
	 */
	public int getScrollValue(ScrollGroup scrollGroup) {
		if (scrollBarsMap.containsKey(scrollGroup)) {			
			ScrollGroupBar bar = scrollBarsMap.get(scrollGroup);
			
			if (bar.isVisible()) {
				return bar.getValue();
			}
		} 			
		return 0;
	}
}
