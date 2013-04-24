package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutComponent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;

/**
 * Single track inside a {@link GBrowserView}. Typically multiple instances
 * are used to construct what user perceives as a track. 
 */
public abstract class Track implements AreaResultListener, LayoutComponent {

	private static final int NAME_VISIBLE_VIEW_RATIO = 20;
	protected GBrowserView view;
	protected List<AreaRequestHandler> areaRequestHandlers = new LinkedList<AreaRequestHandler>();
	protected Strand strand = Strand.FORWARD;
	protected int layoutHeight;
	protected boolean visible = true;
	protected LayoutMode layoutMode = LayoutMode.FIXED;
	protected LayoutMode defaultLayoutMode = LayoutMode.FIXED;	

	public void setView(GBrowserView view) {
    	this.view = view;
    }
    
    /**
     * @param areaRequestHandler
     * @return index of added {@link AreaRequestHandler}
     */
    public int addAreaRequestHandler (AreaRequestHandler areaRequestHandler) {
    	this.areaRequestHandlers.add(areaRequestHandler);
    	return areaRequestHandlers.size() - 1;
    }

	/**
	 * Should be called after Track object is created, but can't be merged to constructor, because the coming areaResult event could cause
	 * call to track object before it's constructed.
	 */
	public void initializeListener() {
		if (areaRequestHandlers != null) {
			for (AreaRequestHandler handler : areaRequestHandlers) {
				view.getQueueManager().addResultListener(handler, this);
			}
		} 
	}

	/**
	 * The method where the actual work of a track typically happens. Each track needs to manage drawables, possibly
	 * caching them.
	 */
	public abstract  Collection<Drawable> getDrawables();
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
	    return areaRequestHandlers != null;
	}
	
    /**
     * Get a map of data sources and column types that this
     * track needs to operate.
     * 
     * Can also return null if this track does not need any data, or column set 
     * can be empty, if the data layer sends the required data anyway (in case of 
     * Conversion classes). 
     */
    public Map<AreaRequestHandler, Set<ColumnType>> requestedData() {
    	
    	return null;
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
	public int getHeight() {
	    return Math.max(layoutHeight, getMinHeight());
	}
	
	/**
	 * Set height of this track.
	 */
    public void setHeight(int height) {
        this.layoutHeight = height;
    }
	
    /**
     * Determine if the track is visible.
     * 
     * @return false.
     */
    public boolean isVisible() {
        return visible;
    }
    
    /**
     * Set track visibility.
     */
    public void setVisible(boolean visible) {
        this.visible = visible;
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
	private String name = "Track";
	private int fullHeight;
	private int FULL_HEIGHT_MARGIN = 10;

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
	
	public String getName() {
		return name ;
	}
	
	public void setName(String name) {
		this.name = name;
	}
	
	public boolean isNameVisible(Rectangle rect) {
		return rect.width > getView().getWidth()/NAME_VISIBLE_VIEW_RATIO;
	}

	protected void drawTextAboveRectangle(String text, Collection<Drawable> drawables, Rectangle rect, int offset) {
		drawables.add(new TextDrawable(rect.x < 0 ? 0 : rect.x, rect.y + offset, text, Color.DARK_GRAY));
	}

	public int getMinHeight() {
		return 0;
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
	}

	public int getFullHeight() {
		if (getLayoutMode() == LayoutMode.FIXED || getLayoutMode() == LayoutMode.FILL) {
			return getHeight();
		} else {
			return Math.max(fullHeight, this.getHeight());
		}
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
}
