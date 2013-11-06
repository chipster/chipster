package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.RenderingHints;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Point2D;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.JComponent;
import javax.swing.JScrollBar;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.SelectionManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

/**
 * Single track inside a {@link GBrowserView}. Typically multiple instances
 * are used to construct a TrackGroup what user perceives as a track. 
 */
public abstract class Track implements DataResultListener, MouseListener {
	
	private JComponent component = new JComponent() {

		@Override
		public void paintComponent(Graphics g) {
			
			//System.out.println("" + this + g);
			
			printTime(null);
			
			//super.paintComponent(g);
			
			Graphics2D g2 = (Graphics2D) g;
			g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,
					RenderingHints.VALUE_ANTIALIAS_ON);

			g2.setPaint(this.getBackground());
			g2.fillRect(0, 0, getWidth(), getHeight());

			// prepare coordinates
			int x = 0;

			List<Selectable> selectables = getSelectables();	
			
			updateSelections(getSelectionManager());
			
 			for (Selectable selectable : selectables) {
				
				List<Drawable> drawables = selectable.getDrawables();
				drawDrawables(g2, x, drawables);
			}
 			
			printTime(getTrackName() + "\tdraw " + selectables.size() + " drawables");
			
		}

		private void drawDrawables(Graphics2D g2, int x,
				Collection<Drawable> drawables) {		
			
			// decide if we will expand drawable for this track
			boolean expandDrawables = canExpandDrawables();

			TrackContext trackContext = null;
			Track thisTrack = Track.this;
			// create view context for this track only if we will use it
			// currently only used for tracks that contain information
			// about reads
			if (expandDrawables && 
					(thisTrack instanceof CoverageTrack ||
							thisTrack instanceof CoverageAverageTrack ||
							thisTrack instanceof CoverageEstimateTrack ||
							thisTrack instanceof QualityCoverageTrack)) {

				if (view.parentPlot.getReadScale() == ReadScale.AUTO) {
					trackContext = new TrackContext(thisTrack, getMaxY(drawables));
				} else {
					trackContext = new TrackContext(thisTrack, view.parentPlot.getReadScale().numReads);
				}
			}
			
			if (LayoutMode.FULL == getLayoutMode()) {
				setFullHeight(drawables);
			}
				
			//Add track height before the track drawables are drawn, because the track coordinates start
			//from the bottom and grow upwards
			int y = convertGraphicsCoordinateToTrack(0);

			for (Drawable drawable : drawables) {

				if(drawable == null) {
					continue;
				}

				// expand drawables to stretch across all height if necessary
				if (expandDrawables) {
					drawable.expand(trackContext);
				}									

				// recalculate position for reversed strands
				int maybeReversedY = (int) y;
				if (isReversed()) {
					maybeReversedY -= getHeight();
				} else {
					drawable.upsideDown();
				}			

				// draw the drawable to the buffer
				drawable.draw(g2, x, maybeReversedY);
			}
		}
		
		@Override
		public Dimension getPreferredSize() {
			Dimension size = super.getPreferredSize();
			if (LayoutMode.FULL == getLayoutMode()) {
				size.height = getFullHeight();	
			} else {
				size.height = getTrackHeight();
			}
			return size;
		}
	};

	private static final int NAME_VISIBLE_VIEW_RATIO = 20;

	private static final float MAX_SELECTION_DISTANCE = 10;
	protected GBrowserView view;
	protected List<DataThread> dataThreads = new LinkedList<DataThread>();
	protected Strand strand = Strand.FORWARD;

	protected LayoutMode layoutMode = LayoutMode.FIXED;
	protected LayoutMode defaultLayoutMode = LayoutMode.FIXED;
	private Map<DataThread, Set<DataType>> dataTypeMap = new HashMap<DataThread, Set<DataType>>();
	private TrackGroup trackGroup;
	
	private long maxViewLength = Long.MAX_VALUE;
	private long minViewLength = 0;
	private int trackHeight = 100;
	
	private String name = "Track";
	private int fullHeight;
	private int FULL_HEIGHT_MARGIN = 10;

	public Track() {
		component.setBackground(Color.white);
		component.setInheritsPopupMenu(true);
		component.addMouseListener(this);
	}

	protected SelectionManager getSelectionManager() {
		return getView().getPlot().getSelectionManager();
	}

	public Track(int height) {
		this();
		this.trackHeight = height;
	}

	public void setView(GBrowserView view) {
    	this.view = view;
		component.addMouseListener(view);
		component.addMouseMotionListener(view);
		component.addMouseWheelListener(view);
    }
    
    /**
     * @param dataThread
     * @return index of added {@link DataThread}
     */
    public int addDataThread (DataThread dataThread) {
    	this.dataThreads.add(dataThread);
    	return dataThreads.size() - 1;
    }

	/**
	 * Should be called after Track object is created, but can't be merged to constructor, because the coming dataResult event could cause
	 * call to track object before it's constructed.
	 */
	public void initializeListener() {
		if (dataThreads != null) {
			for (DataThread handler : dataThreads) {
				view.getQueueManager().addDataResultListener(handler, this);
			}
		} 
	}

	private Integer getMaxY(Collection<Drawable> drawables) {

		int maxY = 0;

		for (Drawable drawable : drawables) {
			maxY = Math.max(drawable.getMaxY() + 1, maxY);            
		}

		return maxY;
	}


	private long time = 0;

	private void printTime(String operation) {
		if (operation == null) {
			time = System.currentTimeMillis();
		} else {
			long elapsed = System.currentTimeMillis() - time;
			time = System.currentTimeMillis();
			if (elapsed > 5) {
				//System.out.println(operation + "\t" + elapsed);
			}
		}
	}

	/**
	 * The method where the actual work of a track typically happens. Each track needs to manage drawables, possibly
	 * caching them.
	 */
	public Collection<Drawable> getDrawables() {
		return null;
	}
	/**
	 * The view under which this track operates.
	 */
	protected GBrowserView getView() {
		return view;
	}
	
	/**
	 * Check if this track has data.
	 */
	public boolean hasData() {
	    return dataThreads != null;
	}
	
    /**
     * Define data sources and dataTypes that this
     * track needs to operate by calling addDataRequest methods.
     * 
     * This empty default implementation is fine if the track doesn't need any data.
     */
    public void defineDataTypes() {
    	
    };
    
    /**
     * Request data of the specified DataType from the first dataThread. Use this method when there is only
     * one dataThread for this track, otherwise use other other overloaded versions to define the
     * dataThreads explicitly.
     * 
     * @param type
     */
    public void addDataType(DataType type) {

    	addDataType(dataThreads.get(0), type);
    }

    /**
     * Request data of the specified DataType from the given dataThread. Use this method when there are several dataThreads 
     * for this track to define the dataThreads explicitly.
     * 
     * @param type
     */
    public void addDataType(DataThread dataThread, DataType type) {

    	Set<DataType> set = dataTypeMap.get(dataThread);
    	
    	if (set == null) {
    		set = new HashSet<DataType>();
    		dataTypeMap.put(dataThread, set);
    	}
    	
    	set.add(type);
    }

	/**
	 * Utility method, return empty Drawable collection.
	 */
	public Collection<Drawable> getEmptyDrawCollection() {
		return new LinkedList<Drawable>();
	}
	
	/**
	 * @return height of this track in pixels.
	 */
	public int getTrackHeight() {
	    return trackHeight;
	}
	
    /**
     * Determine if the track is visible.
     * 
     * @return false.
     */
    public boolean isSuitableViewLength() {
    	
        return (getView().getBpRegion().getLength() >= minViewLength &&
                getView().getBpRegion().getLength() < maxViewLength);
    }   

	public void setStrand(Strand s) {
		this.strand = s;
	}

	public Strand getStrand() {
		return strand;
	}
	
	/**
	 * Determine if this track represents a reverse strand.
	 * 
	 * @return true if this track represents a reverse strand,
	 * false otherwise.
	 */
    public boolean isReversed() {
        return strand == Strand.REVERSE;
    }
    
    /**
     * Determine if drawable elements inside this track can be
     * expanded to stretch across all available height. 
     */
    public boolean canExpandDrawables() {
        return false;
    }

	private Point2D[] arrowPoints = new Point2D[] { new Point.Double(0, 0.25), new Point.Double(0.5, 0.25), new Point.Double(0.5, 0), new Point.Double(1, 0.5), new Point.Double(0.5, 1), new Point.Double(0.5, 0.75), new Point.Double(0, 0.75), new Point.Double(0, 0.25) };

	/**
	 * DOCME
	 * 
	 * Utility method for gettting drawable objects for an arrow figure. The x and y parameters set 
	 * the upper left corner coordinates of the figure and width and height the size of the figure.
	 * This method was placed here to be able to use similar arrows in every track, but more
	 * versatile and separate collection of drawable utilities might be better place in the future.
	 * 
	 * @param x
	 * @param y
	 * @param width
	 * @param height
	 * @return
	 */
	protected Collection<? extends Drawable> getArrowDrawables(int x, int y, int width, int height) {

		Collection<Drawable> parts = getEmptyDrawCollection();

		for (int i = 1; i < arrowPoints.length; i++) {
			Point2D p1 = arrowPoints[i - 1];
			Point2D p2 = arrowPoints[i];

			Point2D p1Scaled = new Point.Double(x + p1.getX() * width, y + p1.getY() * height);
			Point2D p2Scaled = new Point.Double(x + p2.getX() * width, y + p2.getY() * height);

			parts.add(new LineDrawable((int) p1Scaled.getX(), (int) p1Scaled.getY(), (int) p2Scaled.getX(), (int) p2Scaled.getY(), Color.black));
		}

		return parts;
	}
	
	public String getTrackName() {
		return name ;
	}
	
	public void setTrackName(String name) {
		this.name = name;
	}
	
	public boolean isNameVisible(Rectangle rect) {
		boolean wideEnough = rect.width > getView().getWidth()/NAME_VISIBLE_VIEW_RATIO;
		return wideEnough;
	}

	protected void drawTextAboveRectangle(String text, Collection<Drawable> drawables, Rectangle rect, int offset) {
		int x = rect.x < 0 ? 0 : rect.x;
		int y = rect.y + offset;
		drawables.add(new TextDrawable(x, y, text, Color.DARK_GRAY));
	}

	public void setFullHeight(Collection<Drawable> drawables) {
		int maxY = 0;
		for (Drawable drawable : drawables) {
			if (drawable.getMaxY() > maxY && view.getWidth() > drawable.x) {
				maxY = drawable.getMaxY();
			}
		}
		
		if (getLayoutMode() == LayoutMode.FULL) {
			maxY += FULL_HEIGHT_MARGIN;
		}
		this.fullHeight = maxY;
		
		//Repaint until the fullHeigth doesn't change
		if (this.fullHeight != getComponent().getHeight()) {
			getComponent().revalidate();	
		}
	}

	public int getFullHeight() {
		return Math.max(fullHeight, Math.max(component.getHeight(), getTrackHeight()));
	}
	
	public void setLayoutMode(LayoutMode mode) {
		this.layoutMode = mode;
	}
	
	public void setDefaultLayoutMode(LayoutMode mode) {
		this.defaultLayoutMode = mode;
	}
	
	public void setDefaultLayoutMode() {
		this.layoutMode = this.defaultLayoutMode;
	}
	
	public LayoutMode getLayoutMode() {
		return layoutMode;
	}

	public void clearDataTypes() {
		this.dataTypeMap.clear();
	}

	public Map<DataThread, Set<DataType>> getDataTypeMap() {
		return dataTypeMap;
	}
	
	/**
	 * Set limits for the track visibility. 
	 * 
	 * @param minViewLength
	 * @param maxViewLength
	 */
	public void setViewLimits(long minViewLength, long maxViewLength) {
		this.minViewLength = minViewLength;
		this.maxViewLength = maxViewLength;
	}

	public TrackGroup getTrackGroup() {
		return trackGroup;
	}

	public void setTrackGroup(TrackGroup trackGroup) {
		this.trackGroup = trackGroup;
	}

	public void setVisible(boolean state) {
		component.setVisible(state);
	}

	public boolean isVisible() {
		return component.isVisible();
	}

	public JComponent getComponent() {
		return component;
	}

	public List<Selectable> getSelectables() {
		Collection<Drawable> drawables = getDrawables();
		PassiveItem item = new PassiveItem((List<Drawable>) drawables);
		LinkedList<Selectable> itemList = new LinkedList<>();
		itemList.add(item);
		return itemList;
	}

	@Override
	public void mouseClicked(MouseEvent e) {
		final Point point = e.getPoint();
		point.y = convertGraphicsCoordinateToTrack(point.y);
		
		//Collect all selectables inside the selection range
		LinkedList<Selectable> closestSelectables = new LinkedList<>();
		
		for (Selectable selectable : getSelectables()) {

			double distance = selectable.getDistance(point);
			if (distance < MAX_SELECTION_DISTANCE) {
				closestSelectables.add(selectable);
			}
		}
		
		//Sort closest first
		Collections.sort(closestSelectables, new Comparator<Selectable>() {
			@Override
			public int compare(Selectable s1, Selectable s2) {
				double d1 = s1.getDistance(point);
				double d2 = s2.getDistance(point);
				
				if (d1 < d2) {
					return -1;					
				} else if (d1 > d2) {
					return 1;				
				} else {
					return 0;
				}				
			}
		});
		
		
		SelectionManager selectionManager = getSelectionManager();
		
		//Clear selection
		if (!e.isControlDown()) {			
			//If the click did not hit anything, just set selection to null 
			selectionManager.set(getDataUrl(), null);
		}		
				
		//There may be one or more selectables with equal distance. Select them all			
		for (Selectable selectable : closestSelectables) {				
			if (selectable.getDistance(point) > closestSelectables.getFirst().getDistance(point)) {
				//Rest of the selectables are further away
				break;
			}

			//Add one more item to selection
			selectionManager.toggle(getDataUrl(), selectable);				
		}		
		
		updateSelections(selectionManager);	
		getView().redraw();
	}

	private void updateSelections(SelectionManager selectionManager) {
		for (Selectable selectable : getSelectables()) {
			selectable.setSelected(selectionManager.isSelected(getDataUrl(), selectable));
		}
	}


	/**
	 * Track coordinates grow upwards, whereas Java Graphics coordinates grow downwards.
	 * 
	 * @param y
	 * @return
	 */
	private int convertGraphicsCoordinateToTrack(int y) {
		
		int height = getComponent().getHeight();
		
		return height - y;
	}

	@Override
	public void mouseEntered(MouseEvent e) {
	}


	@Override
	public void mouseExited(MouseEvent e) {
	}
	

	@Override
	public void mousePressed(MouseEvent e) {
	}


	@Override
	public void mouseReleased(MouseEvent e) {
	}

	public DataUrl getDataUrl() {
		if (dataThreads.size() > 0) {
			return dataThreads.get(0).getDataSource().getDataUrl();
		}
		return null;
	}

	public boolean isShowMoreCapable() {
		return false;
	}

	public int getVisibleWidth() {
		JScrollBar bar = getTrackGroup().getScrollGroup().getComponent().getVerticalScrollBar();
		int width = view.getWidth();
	
		if (bar.isVisible()) {
			return width - bar.getWidth();
		} else {
			return width;
		}
	}
}
