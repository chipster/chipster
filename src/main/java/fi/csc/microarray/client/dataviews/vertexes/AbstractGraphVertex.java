package fi.csc.microarray.client.dataviews.vertexes;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Point;
import java.awt.Stroke;
import java.awt.geom.Rectangle2D;
import java.util.Hashtable;
import java.util.Map;

import javax.swing.BorderFactory;
import javax.swing.border.Border;

import org.jgraph.graph.DefaultGraphCell;
import org.jgraph.graph.GraphConstants;

import fi.csc.microarray.client.dataview.MicroarrayGraph;
import fi.csc.microarray.client.dataview.PositionChangeAnimator;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataItem;

/**
 * Abstract class for vertexes used on MicroarrayGraph
 * 
 * @author mkoski
 * 
 */
public abstract class AbstractGraphVertex extends DefaultGraphCell {

	private MicroarrayGraph graph;
	
	public static Color DEFAULT_VERTEX_COLOR = VisualConstants.CATEGORY_COLORS[0]; 
	public static final Border BORDER_NORMAL =
		BorderFactory.createLineBorder(Color.black);
	
	public static final Stroke SELECTION_STROKE = new BasicStroke(2);
	public static final Color SELECTED_BORDER_COLOR = Color.BLACK;
	
	public static final int DEFAULT_WIDTH = 40;
	public static final int DEFAULT_HEIGHT = 20;
	
	public static final int DEFAULT_MARGIN_X = 10;
	public static final int DEFAULT_MARGIN_Y = 30; 

	/**
	 * Constructor for graph
	 * 
	 * @param x
	 * @param y
	 * @param data
	 * @param graph
	 */
	public AbstractGraphVertex(int x, int y, DataItem data,
			MicroarrayGraph graph) {
		super(data);
		this.graph = graph;
		this.setPosition(new Point(x,y));
		
		GraphConstants.setGradientColor(
				this.getAttributes(),
				DEFAULT_VERTEX_COLOR);
		
		GraphConstants.setBorderColor(
				this.getAttributes(), 
				Color.black);
		
		GraphConstants.setOpaque(
				this.getAttributes(), 
				true);
	}

	/**
	 * horizontal position
	 * 
	 * @return
	 */
	public int getX() {
		return (int) this.getBounds().getX();
	}

	/**
	 * vertical position
	 * 
	 * @return
	 */
	public int getY() {
		return (int) this.getBounds().getY();
	}

	/**
	 * Sets new value to horizontal and vertical positions
	 * @param x
	 * @param y
	 */
	public void setPosition(Point point) {

		Map attrs = new Hashtable();
		GraphConstants.setBounds(attrs, new Rectangle2D.Double(
				point.getX(), point.getY(), getDefaultWidth(), getDefaultHeight()));
				
		graph.getGraphLayoutCache().editCell(this, attrs);
		
		graph.repaint();
	}
	
	public void setPositionAnimated(Point target) {
		new PositionChangeAnimator(
				this,								 //what
				new Point(this.getX(), this.getY()), //from
				target);							 //to		
	}

	/**
	 * Gets the vertex bounds
	 * @return
	 */
	public Rectangle2D getBounds() {
		return GraphConstants.getBounds(this.getAttributes());
	}
	
	/**
	 * Gets the graph object
	 * @return
	 */
	public MicroarrayGraph getGraph() {
		return this.graph;
	}

	/**
	 * Gets default width of the vertex
	 * @return
	 */
	public abstract int getDefaultWidth();

	/**
	 * Gets default height of the vertex
	 * @return
	 */
	public abstract int getDefaultHeight();

	/**
	 * Overrides toString method. The output of this method is used
	 * as a label of the vertex
	 */
	@Override
	public abstract String toString();

	/**
	 * Returns tooltip string
	 * @return
	 */
	public abstract String getToolTipString();
	
	public abstract void addChildVertex(GraphVertex child);
	
	public abstract boolean removeChildVertex(GraphVertex child);
	
	
}
