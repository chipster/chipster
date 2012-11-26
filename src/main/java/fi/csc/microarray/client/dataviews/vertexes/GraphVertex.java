package fi.csc.microarray.client.dataviews.vertexes;

import java.awt.Color;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;
import org.jgraph.graph.GraphConstants;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataview.MicroarrayGraph;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * This class defines a vertex of the graph view. Each vertex represents one
 * dataset, and they are hierarchically connected to show the workflows that
 * were needed to produce each set of data. The color of the vertex shows
 * which type of data is in question (that is, what kind of operation steps
 * it has already undergone). The vertices can be selected (which will select 
 * the corresponding dataset) and moved around on the graph panel.
 * 
 * @author Janne Käki, Mikko Koski
 *
 */
public class GraphVertex extends AbstractGraphVertex {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(GraphVertex.class);
	
	private static final String AT_LEAST_ROWS_CACHENAME = "at-least-rows";
	private static final int MAX_ROWS_TO_COUNT = 1000; 
	private static final int MAX_BYTES_TO_COUNT = 100*1024;
	
	private List<GraphVertex> children;
	
	private String tooltipText;
	
	private MicroarrayGraph graph;
	
	/**
	 * Creates a new graph vertex. The coordinates of the new vertex are
	 * defined right and down from the upper left corner of the graph area.
	 * The parent of the vertex is defined by the parent dataset of the
	 * given data parameter.
	 * 
	 * @param x Horizontal position of the new vertex (in pixel coordinates).
	 * @param y Vertical position of the new vertex (in pixel coordinates).
	 * @param data Dataset which this vertex is to represent on the graph.
	 * @throws NullPointerException If the provided DataBean is null.
	 */
	public GraphVertex(int x, int y, DataBean data, MicroarrayGraph graph) {
		super(x, y, data, graph);
		
		children = new ArrayList<GraphVertex>();
		
		Color c;
		if (data == null || data.getOperationRecord() == null || data.getOperationRecord().getCategoryColor() == null) {
				c = DEFAULT_VERTEX_COLOR;			
				
		} else {			
			c = data.getOperationRecord().getCategoryColor();
		}
		
		GraphConstants.setGradientColor(this.getAttributes(),c);
		
		GraphConstants.setBackground(this.getAttributes(),c.brighter());
                
		GraphConstants.setOpaque(this.getAttributes(), true);
		GraphConstants.setBorder(this.getAttributes(), BORDER_NORMAL);
		this.addPort();
		
		this.graph = graph;
		
		try {
			int rowCount;

			if (data.getContentLength() < MAX_BYTES_TO_COUNT) {
				// check if rows are counted already 
				Integer cachedCount = (Integer)data.getFromContentBoundCache(AT_LEAST_ROWS_CACHENAME);
				if (cachedCount != null) {
					rowCount = cachedCount; 

				} else {
					// count rows
					Table rowCounter = data.queryFeatures("/column/*").asTable();
					rowCount = 0;
					while (rowCounter != null && rowCounter.nextRow() && rowCount < MAX_ROWS_TO_COUNT) {
						rowCount++;
					}
					data.putToContentBoundCache(AT_LEAST_ROWS_CACHENAME, (Integer)rowCount);
				}

				if (rowCount < 1000) {
					tooltipText = ", " + rowCount  + " rows";
				} else {
					tooltipText = ", over " + rowCount  + " rows";
				}
				
			} else {
				tooltipText = ""; // too large file, maybe binary, cannot say anything 
			}
			
		} catch (MicroarrayException e) {
			tooltipText = "";
		}
		
		graph.getGraphLayoutCache().toFront(new Object[] { this } );
		
		graph.repaint();
		
	}
	
	/**
	 * Gets the user specified data for this vertex. Data is stored in the 
	 * construction method
	 * 
	 * @return user specified data
	 */
	public DataBean getData() {
		Object userObject = this.getUserObject();
		assert (userObject instanceof DataBean);
		return (DataBean) userObject;
	}

	/**
	 * Gets the folder where the data belongs to
	 * @return
	 */
	public DataFolder getFolder() {
		return this.getData().getParent();
	}
	
	/**
	 * Gets the group of the vertex. If vertex does not have a group (is not raw data)
	 * returns null.
	 *
	 * @return group or <code>null</code> if vertex is not a member of any group
	 */
	public GroupVertex getGroup(){
		if (this.getParent() instanceof GroupVertex) {
			return (GroupVertex)this.getParent();
		} else {
			return null;
		}
	}
	
	public boolean isInGroup(){
		return getGroup() != null;
	}
	
	@Override
	public int getDefaultHeight() {
		return DEFAULT_HEIGHT;
	}

	@Override
	public int getDefaultWidth() {
		return DEFAULT_WIDTH;
	}
	
	public int getMarginX(){
		return DEFAULT_MARGIN_X;
	}
	
	public int getMarginY(){
		return DEFAULT_MARGIN_Y;
	}
	
	/**
	 * Gets the tooltip text of the selected cell
	 * @return tooltip text
	 */
	public String getToolTipString() {
		return this.getData().toString() + tooltipText;
	}
	
	
	/** 
	 * Returns first four letters of the category name or all if the name is shorter
	 * 
	 * @see fi.csc.microarray.client.dataviews.vertexes.AbstractGraphVertex#toString()
	 */
	@Override
	public String toString() {
		return Session.getSession().getPrimaryModule().getShortDataName(getData().getName());		
	}
	
	public void addChildVertex(GraphVertex child){
		children.add(child);
	}
	
	public boolean removeChildVertex(GraphVertex child){
		if(children.contains(child)){
			children.remove(child);
			return true;
		} else {
			return false;
		}
	}
	
	/**
	 * Gets child vertices of the current vertex
	 * @return
	 */
	public List<GraphVertex> getChildVertices() {
		return this.children;
	}
	
	/**
	 * Gets root vertexes of this vertex
	 * 
	 * @return
	 */
	public List<GraphVertex> getRootVerticesOfThisVertex(){
		List<GraphVertex> roots = new ArrayList<GraphVertex>();		
		findRoots(this.getData(), roots);
		return roots;
	}
	
	private void findRoots(DataBean data, List<GraphVertex> roots) {
		for (DataBean source : data.getLinkTargets(Link.DERIVATION, Link.MODIFICATION)) {
			if(source.getLinkTargets(Link.DERIVATION, Link.MODIFICATION).size() == 0) {
				GraphVertex vertex = (GraphVertex)(graph.getVertexMap().get(source));
				if(!roots.contains(vertex)){
					roots.add((GraphVertex)(graph.getVertexMap().get(source)));
				}
			} else {
				this.findRoots(source, roots);
			}
		}
	}
	
	public boolean isRoot() {
		logger.debug("vertex " + this.getData().getName() + " derives from " + this.getData().getLinkTargets(Link.DERIVATION).size() + " databeans");
		return this.getData().getLinkTargets(Link.derivationalTypes()).size() == 0;
	}
}	
