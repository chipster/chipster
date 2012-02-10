package fi.csc.microarray.client.dataview;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.dataviews.vertexes.AbstractGraphVertex;
import fi.csc.microarray.client.dataviews.vertexes.GraphVertex;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;

/**
 * Class for some  methods which are iterating through
 * the vertexes to make placement of the new vertexes as good as possible
 * 
 * @author Petri Klemelä, Aleksi Kallio
 */
public class LayoutManager {
	
	private MicroarrayGraph graph;
	
	private static final Logger logger = Logger.getLogger(LayoutManager.class);
	
	public LayoutManager(MicroarrayGraph graph) {
		this.graph = graph;
	}
	
	/** 
	 * Uses helper methods of this class to find a place for a new vertex and to make
	 * space if it is necessary.
	 * 
	 * @param inserted GraphVertex, which position should be updated
	 */
	public void updateLayout(AbstractGraphVertex inserted) {		

		if (!(inserted instanceof GraphVertex)) {
			throw new IllegalArgumentException("vertex of type " + inserted.getClass().getSimpleName() + " not supported");
		}
		
		GraphVertex vertex = (GraphVertex)inserted;
		
		logger.debug(vertex.getData().getName() + " is root: " + vertex.isRoot());
		
		if (vertex.isRoot()){
			vertex.setPosition(getNewRootPosition());
			
		} else {
			Point preferredPlace = getPreferredPlace(vertex);
			vertex.setPosition(preferredPlace);
		}
		
		logger.debug(vertex.getData().getName() + " got bounds: " + vertex.getBounds());

	}
	
	/** 
	 * Returns the farthest vertex position added by default margins.
	 * Relies on the getGraphSize()-method. 
	 * 
	 * @return the coordinates for the next root
	 */
	public Point getNewRootPosition(){
		Dimension dim = graph.getGraphSize();
		int x = (int)dim.getWidth();
		
		//Snap to left border if close enough
		if(x < 60){
			x = 10;
		}
		return new Point(x, 10);
	}
	
	/** Finds a place for the new cell according it's sources and siblings.
	 * 
	 * @param vertex - cell which is inserted
	 * @return Point where the cell should be placed
	 */
	private Point getPreferredPlace(GraphVertex vertex) {
		
		int leftSourceX = Integer.MAX_VALUE;
		int bottomSourceY = 0;
		
		List<Integer> sourceXs = new ArrayList<Integer>();	

		// check sources
		for (DataBean sourceBean : vertex.getData().getLinkTargets(Link.derivationalTypes())) {										
			GraphVertex source = graph.getVertexMap().get(sourceBean);
			
			// find the farthest source from the left					
			if (source.getX() < leftSourceX) {
				leftSourceX = source.getX();
			}
			
			sourceXs.add(source.getX());

			// find the lowest source
			if (source.getY() > bottomSourceY) {
				bottomSourceY = source.getY();
			}			
		}
				
		// if there is free space below the left source
		int sourceY = bottomSourceY+vertex.getMarginY()+vertex.getDefaultHeight();
		Point potentialPosition = new Point(leftSourceX,sourceY);

		if (getIntersectingVertex(potentialPosition, vertex) == null) {
			return potentialPosition;
		}
		
		// or if below any other source
		for (int sourceX : sourceXs) {
			Point point = new Point(sourceX, sourceY);
			if (getIntersectingVertex(point, vertex) == null) {
				return point;
			}
		}
		
		AbstractGraphVertex intersectingCell;
		boolean done = false;
		
		// start below left source and go to the right until the free place is found
		while (!done) {			
			intersectingCell = this.getIntersectingVertex(potentialPosition, vertex);
			
			if (intersectingCell != null && intersectingCell instanceof GraphVertex) {
				potentialPosition = new Point((int)(potentialPosition.getX() + vertex.getDefaultWidth() + vertex.getMarginX()), (int)(potentialPosition.getY()));
				
			} else {
				done = true;
			}						
		}
		return potentialPosition;
	}
	
	private AbstractGraphVertex getIntersectingVertex(Point newPoint, GraphVertex vertex) {
		
		Dimension dim = new Dimension(GraphVertex.DEFAULT_WIDTH,GraphVertex.DEFAULT_HEIGHT);
		Rectangle newPlace = new Rectangle(newPoint, dim);
		
		for (AbstractGraphVertex cell : graph.getVisibleVertexes()) {
			logger.debug("does " + cell + " intersect " + vertex.getData().getName() + ": " + (cell.getBounds().intersects(newPlace)  && cell != vertex));  
			if (cell.getBounds().intersects(newPlace) && cell != vertex) {
				return cell;
			}
		}
		return null;
	}	
}
