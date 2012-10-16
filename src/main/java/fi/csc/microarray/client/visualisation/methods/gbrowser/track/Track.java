package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.geom.Point2D;
import java.util.Collection;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;

/**
 * Single track inside a {@link View}. Typically multiple instances
 * are used to construct what user perceives as a track. 
 */
public abstract class Track implements AreaResultListener {

	private static final int NAME_VISIBLE_VIEW_RATIO = 20;
	protected View view;
	protected DataSource file;
	protected Strand strand = Strand.FORWARD;
	protected Integer height;
	protected boolean visible = true;
	
    public Track(View view, DataSource file) {
		this.view = view;
		this.file = file;
	}

	/**
	 * Should be called after Track object is created, but can't be merged to constructor, because the coming areaResult event could cause
	 * call to track object before it's constructed.
	 */
	public void initializeListener() {
		if (file != null) {
			view.getQueueManager().addResultListener(file, this);
		} else {
			throw new RuntimeException("Track has no file: " + this);
		}
	}

	/**
	 * The method where the actual work of a track typically happens. Each track needs to manage drawables, possibly
	 * caching them.
	 */
	public abstract Collection<Drawable> getDrawables();

	/**
	 * The view under which this track operates.
	 */
	protected View getView() {
		return view;
	}
	
	/**
	 * Check if this track has data.
	 */
	public boolean hasData() {
	    return file != null;
	}
	
    /**
     * Get a map of data sources and column types that this
     * track needs to operate.
     * 
     * Can also return null if this track does not need any data.
     */
	public abstract Map<DataSource, Set<ColumnType>> requestedData();

	/**
	 * If track is concised, it is not showing exact data but approximations calculated from the data.
	 */
	public abstract boolean isConcised();

	/**
	 * Utility method, return empty Drawable collection.
	 */
	public Collection<Drawable> getEmptyDrawCollection() {
		return new LinkedList<Drawable>();
	}
	
	/**
	 * Each track has individual height. If it is not set explicitly,
	 * the default height is taken from View.
	 * 
	 * @return height of this track in pixels.
	 */
	public Integer getHeight() {
	    return height;
	}
	
	/**
	 * Set height of this track.
	 */
    public void setHeight(Integer height) {
        this.height = height;
    }

	/**
	 * Determine if the track can be resized vertically.
	 * 
	 * @return true if track can be resized, false if it has
	 * static height.
	 */
	public abstract boolean isStretchable();
	
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
        return strand == Strand.REVERSED;
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
}
