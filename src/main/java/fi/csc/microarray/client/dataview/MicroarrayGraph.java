package fi.csc.microarray.client.dataview;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Point;
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Vector;

import javax.swing.ToolTipManager;

import org.apache.log4j.Logger;
import org.jgraph.JGraph;
import org.jgraph.event.GraphModelEvent;
import org.jgraph.event.GraphModelListener;
import org.jgraph.graph.CellView;
import org.jgraph.graph.DefaultEdge;
import org.jgraph.graph.DefaultPort;
import org.jgraph.graph.GraphConstants;
import org.jgraph.graph.GraphLayoutCache;
import org.jgraph.graph.GraphModel;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataviews.vertexes.AbstractGraphVertex;
import fi.csc.microarray.client.dataviews.vertexes.GraphVertex;
import fi.csc.microarray.client.dataviews.vertexes.GroupVertex;
import fi.csc.microarray.client.dataviews.vertexes.PhenodataVertex;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.databeans.ContentChangedEvent;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataItemCreatedEvent;
import fi.csc.microarray.databeans.DataItemRemovedEvent;
import fi.csc.microarray.databeans.LinksChangedEvent;
import fi.csc.microarray.util.SwingTools;

/**
 * A graph element for workflow view based on JGraph.
 * 
 * @author Mikko Koski, Aleksi Kallio
 * 
 */
public class MicroarrayGraph extends JGraph implements DataChangeListener, PropertyChangeListener, GraphModelListener {

	private static final Logger logger = Logger.getLogger(MicroarrayGraph.class);

	private ClientApplication application = Session.getSession().getApplication();

	/**
	 * Maps dataset vertexes and datas together
	 */
	private Map<DataBean, GraphVertex> vertexMap = new HashMap<DataBean, GraphVertex>();
	private List<GroupVertex> groups = new ArrayList<GroupVertex>();
	private LayoutManager layoutManager = new LayoutManager(this);
	private GraphPanel graphPanel;
	private GraphModel model;

	public MicroarrayGraph(GraphModel model, GraphLayoutCache cache, GraphPanel graphPanel) {
		super(model, cache);

		this.model = model;
		this.graphPanel = graphPanel;
		this.setBackground(Color.WHITE);

		// Direct call to parents L&F update, see this.updateUI()
		super.updateUI();

		/*
		 * Regardles of the updateUI-call jus before the font size of graph cell is fixed only when the first cell is added. To handle
		 * situations when the font size is changed before importing any data we force the jgraph to initialise itself by inserting dummy
		 * vertex and removing it. This *should* happen before any repaint operation, and thus the vertex won't ever show.
		 */
		GroupVertex vertex = new GroupVertex(0, 0, null, this);
		getGraphLayoutCache().insert(vertex);
		getGraphLayoutCache().remove(new Object[] { vertex });

		// registers the graph object with the tooltip manager
		ToolTipManager.sharedInstance().registerComponent(this);
		
		model.addGraphModelListener(this);

		// start listening
		application.getDataManager().addDataChangeListener(this);
		application.addClientEventListener(this);
	}

	/**
	 * Inserts a new DataBean to this graph.
	 * 
	 * @param data
	 *            New DataBean to be inserted to this view.
	 */
	public void insertData(DataBean data) {

		// check parameters
		if (vertexMap.containsKey(data)) {
			throw new RuntimeException(data.getName() + " was already present in the graph");
		}
		
		boolean savedPosition = data.getX() != null && data.getY() != null; 
		
		int x = 0;
		int y = 0;
		
		if (savedPosition) {
			x = data.getX();
			y = data.getY();
		}		

		// create approriate type of vertex
		GraphVertex vertex;
		if (data.queryFeatures("/phenodata").exists()) {
			// create a metadata vertex
			vertex = new PhenodataVertex(x, y, data, this);			

		} else if (data.getLinkTargets(Link.derivationalTypes()).size() == 0) {
			// create a root vertex
			vertex = new GraphVertex(x, y, data, this);		

		} else {
			// create a normal vertex
			if (savedPosition) {
				vertex = new GraphVertex(x, y, data, this);
			} else {
				vertex = new GraphVertex(10, 10, data, this);
			}
		}
		
		vertex.setAllowsAutoLayout(!savedPosition);

		// create derivational (DERIVATION/MODIFICATION) links
		for (Link type : Link.derivationalTypes()) {
			for (DataBean targetBean : data.getLinkTargets(type)) {
				GraphVertex parent = vertexMap.get(targetBean);
				// Note that PARENT is the TARGET of derivational link and the
				// child is the SOURCE
				insertLink(vertex, parent, type, data);
			}
		}

		// make vertex visible
		vertexMap.put(data, vertex);
		getGraphLayoutCache().insert(vertex);

		// The placement is unknown before the links are created. These stay valid only
		// if no derivation/modification links are added.
		
		layoutManager.updateLayout(vertex, data);
		
		graphPanel.autoZoom();
		scrollCellToVisibleAnimated(vertex);
		repaint();
	}

	/**
	 * Removes the given dataset from this graph.
	 * 
	 * @param data
	 *            Dataset to be removed.
	 * @return True if remove succeeded, false if it failed (that is, corresponding vertex was not found in the graph).
	 */
	public void removeData(DataBean data) {

		// fetch vertex
		GraphVertex vertex = vertexMap.get(data);
		if (vertex == null) {
			throw new IllegalArgumentException(data.getName() + " was not present in graph");
		}

		// mark vertex to be removed
		Vector<Object> removedCells = new Vector<Object>();
		removedCells.add(vertex);

		// process children and their nested objects to be removed (edges and ports)
		List children = vertex.getChildren();
		for (Object child : children) {
			if (child instanceof DefaultPort) {
				DefaultPort port = (DefaultPort) child;
				Set edges = port.getEdges();

				for (Object edgeObject : edges) {
					if (edgeObject instanceof DefaultEdge) {
						DefaultEdge edge = (DefaultEdge) edgeObject;
						Object source = edge.getSource();
						if (source instanceof DefaultPort) {
							DefaultPort sourcePort = (DefaultPort) source;
							removedCells.add(sourcePort);
						}
						removedCells.add(edge);
					}
				}

				removedCells.add(port);
			}
		}

		// removes the group vertex if there is no group members left in the group
		if (vertex.getGroup() != null) {
			GroupVertex group = vertex.getGroup();

			// if there are no other vertexes in the group remove it
			if ((group.getChildCount() - 2) < 1) {
				removedCells.add(group);
			}
		}

		// remove from graph and graph model
		this.getGraphLayoutCache().remove(removedCells.toArray()); // removed from visible graph (?)
		this.model.remove(removedCells.toArray()); // removes from graph model

		// remove from vertex map
		vertexMap.remove(data);

		// remove from other structures
		for (Object o : removedCells) {
			if (groups.contains(o)) {
				// Data is a group
				groups.remove(o); // remove from groups list

			} else if (o instanceof AbstractGraphVertex) {
				// Data is not a root and not a group => remove the data from parent's list of children
				Object parent = ((AbstractGraphVertex) o).getParent();
				if (parent != null) {
					if (parent instanceof GraphVertex) {
						((GraphVertex) parent).removeChildVertex((GraphVertex) o);

					} else {
						// GroupVertex does not store information about it's
						// children. This information is stored only by JGraph
						// and it is removed when cell is removed from
						// graphLayoutCache. So we do nothing
					}
				}
			}
		}

		graphPanel.autoZoom();
		this.repaint();
	}

	public void createGroup(DataBean data) {
		GraphVertex groupMember = vertexMap.get(data);
		createGroup(new GraphVertex[] { groupMember });
	}

	private GroupVertex createGroup(GraphVertex[] children) {
		if (children.length < 1) {
			throw new IllegalArgumentException("vertex list is empty");
		}

		// Gets children bounds so group cell can be positioned to upper left corner
		Rectangle2D bounds = this.getCellBounds(children);
		Point p = new Point((int) bounds.getX(), (int) bounds.getY());
		// Point p = layoutManager.getNewRootPosition();

		GroupVertex group = new GroupVertex((int) p.getX(), (int) p.getY(), children, this);
		this.getGraphLayoutCache().insertGroup(group, children);

		group.collapse();

		groups.add(group);

		return group;
	}

	public void scrollCellToVisibleAnimated(GraphVertex cell) {
		if (cell != null) {
			new ScrollAnimator(graphPanel, cell.getBounds().getBounds());
		}
	}

	/**
	 * Updates graph's selected vertexes. This is done by synchronizing selected datas and graph cells. Update should be done when
	 * selections are changes (of course) and also if groups are collapsed or expanded.
	 * 
	 */
	public void updateSelectedCells() {
		graphPanel.setInternalSelection(true);

		List<Object> selectedCells = new ArrayList<Object>();

		for (DataBean selected : application.getSelectionManager().getSelectedDatasAsArray()) {
			selectedCells.add(this.getVertexMap().get(selected));
		}

		this.setSelectionCells(selectedCells.toArray());

		// If group is collapsed and all of its content is selected the group should
		// be selected too
		for (GroupVertex group : groups) {
			if (group.isAllChildrenSelected()) {
				this.addSelectionCell(group);

			} else {
				this.removeSelectionCell(group);
			}
		}

		this.repaint();

		graphPanel.setInternalSelection(false);
	}

	/**
	 * Returns the tooltip to be displayed for cell gestured by the event. If vertex was not gestured, returns null.
	 */
	@Override
	public String getToolTipText(MouseEvent event) {
		Object cell = getFirstCellForLocation(event.getX(), event.getY());
		if (cell instanceof GraphVertex) {
			return ((GraphVertex) cell).getToolTipString();
		}
		return null;
	}

	public Map<DataBean, GraphVertex> getVertexMap() {
		return vertexMap;
	}

	/**
	 * Returns all vertexes, including invisible ones.
	 */
	public List<AbstractGraphVertex> getAllVertexes() {
		List<AbstractGraphVertex> vertexes = new ArrayList<AbstractGraphVertex>();
		for (DataBean data : application.getDataManager().databeans()) {
			GraphVertex vertex = vertexMap.get(data);
			if (vertex != null) {
				vertexes.add(vertex);
			}
		}
		return vertexes;
	}

	/**
	 * Gets all visible vertexes on the graph.
	 */
	public List<AbstractGraphVertex> getVisibleVertexes() {
		CellView[] views = this.getGraphLayoutCache().getCellViews();
		List<AbstractGraphVertex> vertexes = new ArrayList<AbstractGraphVertex>();
		for (CellView view : views) {
			if (view != null && view.getCell() instanceof AbstractGraphVertex) {
				AbstractGraphVertex vertex = (AbstractGraphVertex) view.getCell();
				vertexes.add(vertex);
			}
		}
		return vertexes;
	}

	/**
	 * Returns vertexes at selected point.
	 */
	public List<AbstractGraphVertex> getVertexesAtPoint(Point2D point) {
		int x = (int) point.getX();
		int y = (int) point.getY();
		List<AbstractGraphVertex> vertexes = new ArrayList<AbstractGraphVertex>();

		// Gets the topmost cell
		AbstractGraphVertex topmost;
		Object current;
		Object temporaryTopmost = getFirstCellForLocation(x, y);

		if (temporaryTopmost instanceof AbstractGraphVertex) {
			topmost = (AbstractGraphVertex) temporaryTopmost;
		} else {
			// If the topmost is not AbstractGraphVertex something is wrong
			return null;
		}

		vertexes.add(topmost);
		current = topmost;

		// adds all vertexes to list while current is not the same as topmost
		do {
			current = getNextCellForLocation(current, x, y);

			if (current != topmost && current instanceof AbstractGraphVertex) {
				vertexes.add((AbstractGraphVertex) current);
			}
		} while (current != topmost);

		return vertexes;
	}

	/**
	 * Gives the dimension from the origo to the farthest vertexes in right and bottom added width default margins of vertexes. This is
	 * different from actual size of the graph, this only tells where the vertexes are situated.
	 */
	public Dimension getGraphSize() {
		Dimension dim = new Dimension(0, 0);

		// find maximum dimensions of vertexes
		for (AbstractGraphVertex vertex : this.getAllVertexes()) {
			int x = vertex.getX();
			int y = vertex.getY();

			if (x > dim.getWidth()) {
				dim.width = x;
			}

			if (y > dim.getHeight()) {
				dim.height = y;
			}
		}

		dim.width += GraphVertex.DEFAULT_WIDTH + GraphVertex.DEFAULT_MARGIN_X;
		dim.height += GraphVertex.DEFAULT_HEIGHT + GraphVertex.DEFAULT_MARGIN_Y;

		return dim;
	}

	/**
	 * Uses JGraph's fromScreen method to convert dimensions from the screen coordinates.
	 */
	public Dimension fromScreenCoordinates(Dimension dim) {
		Point2D fromScreen = this.fromScreen(new Point.Double(dim.getWidth(), dim.getHeight()));
		return new Dimension((int) fromScreen.getX(), (int) fromScreen.getY());
	}

	/**
	 * Gets first cell for location and otherwise that the method of the super class does, this method ignores expanded group cells. If user
	 * clicked vertex which is in a expanded groups bounds the super class method returned cell of the expanded group, which is in fact
	 * useless, because user wants to select the dataset cell, not the group cell.
	 */
	@Override
	public Object getFirstCellForLocation(double x, double y) {
		// get cell
		Object cell = super.getFirstCellForLocation(x, y);

		// Point p = new Point((int)x, (int)y);

		// accept only collapsed groups and graph vertexes
		AbstractGraphVertex vertex;
		if (cell instanceof GroupVertex) {
			// Collapsed group
			vertex = (GroupVertex) cell;
		} else if (cell instanceof GraphVertex) {
			// Graph vertex. OK!
			vertex = (GraphVertex) cell;
		} else {
			vertex = null;
		}
		return vertex;
	}

	/**
	 * Sets UI for graph. The BasicGraphUI does autoscrolling in a really annoying way and it can't be disabled very easily. That's why we
	 * use the own SmartAutoscrollGraphUI
	 * 
	 * @Override public void updateUI() { setUI(new SmartAutoscrollGraphUI()); invalidate(); }
	 */

	public void propertyChange(PropertyChangeEvent event) {
		if (event instanceof DatasetChoiceEvent && event.getSource() != this) {
			SwingTools.runInEventDispatchThread(new Runnable() {
				public void run() {
					updateSelectedCells();		
				}
			});
		}
	}

	/**
	 * Iterates through every vertex and tells them to update their look to show the selected and non selected items
	 */
	public void dataChanged(final DataChangeEvent event) {

		SwingTools.runInEventDispatchThread(new Runnable() {
			public void run() {

				if (event instanceof ContentChangedEvent) {
					// repaint graph if phenodata changed. This adds/removes the warning icon
					// Jgraph 5.12 is very stingy with it's repaints and all following three calls
					// are needed to force repaint
					setVisible(false);
					setVisible(true);
					addOffscreenDirty(getBounds());

				} else if (event instanceof LinksChangedEvent) {
					LinksChangedEvent lcEvent = (LinksChangedEvent) event;
					if (lcEvent.isCreation()) {
						insertLink(lcEvent.getSource(), lcEvent.getTarget(), lcEvent.getType());
					} else {
						removeLink(lcEvent.getSource(), lcEvent.getTarget(), lcEvent.getType());
					}

				} else if (event instanceof DataChangeEvent) {

					DataItem data = ((DataChangeEvent) event).getDataItem();
					if (data instanceof DataBean) {
						// only beans are shown in graph
						if (event instanceof DataItemCreatedEvent) {
							insertData((DataBean) data);

						} else if (event instanceof DataItemRemovedEvent) {
							removeData((DataBean) data);
						}

						repaint();
					}
				}
			}
		});

	}

	private void moveCloseToAnnotated(DataBean data) {
		List<DataBean> list = data.getLinkTargets(DataBean.Link.ANNOTATION);
		DataBean annotatedBean = list.size() > 0 ? list.get(0) : null;
		if (annotatedBean == null) {
			throw new IllegalArgumentException(data.getName() + " must have ANNOTATION link target(s)");
		}

		GraphVertex vertex = vertexMap.get(data);
		Rectangle2D bounds = vertexMap.get(annotatedBean).getBounds();
		vertex.setPosition(new Point((int) bounds.getX() + GraphVertex.DEFAULT_WIDTH + GraphVertex.DEFAULT_MARGIN_X, (int) bounds.getY()));
	}

	private static void setEdgeStyle(DefaultEdge edge) {
		GraphConstants.setSelectable(edge.getAttributes(), false);
	}

	private static void setModificationEdgeStyle(DefaultEdge edge) {
		setEdgeStyle(edge);
		GraphConstants.setLineBegin(edge.getAttributes(), GraphConstants.ARROW_TECHNICAL);
		GraphConstants.setBeginFill(edge.getAttributes(), true);
		GraphConstants.setBeginSize(edge.getAttributes(), 5);
		GraphConstants.setForeground(edge.getAttributes(), Color.LIGHT_GRAY);
		GraphConstants.setDashPattern(edge.getAttributes(), new float[] { 1f, 2f });
	}

	private static void setDerivationEdgeStyle(DefaultEdge edge) {
		setEdgeStyle(edge);
		GraphConstants.setLineBegin(edge.getAttributes(), GraphConstants.ARROW_TECHNICAL);
		GraphConstants.setBeginFill(edge.getAttributes(), true);
		GraphConstants.setBeginSize(edge.getAttributes(), 5);
		GraphConstants.setLineWidth(edge.getAttributes(), 1);
	}

	private static void setAnnotationEdgeStyle(DefaultEdge edge) {
		setEdgeStyle(edge);
		GraphConstants.setForeground(edge.getAttributes(), Color.LIGHT_GRAY);
		GraphConstants.setDashPattern(edge.getAttributes(), new float[] { 2f, 2f });
	}

	/**
	 * Creates a link edge between beans.
	 */
	public void insertLink(DataBean source, DataBean target, Link type) {
		logger.debug("adding link (" + type + ") from " + source.getName() + " to " + target.getName());

		// check vertices
		GraphVertex sourceVertex = vertexMap.get(source);
		GraphVertex targetVertex = vertexMap.get(target);
		if (sourceVertex == null || targetVertex == null) {
			throw new RuntimeException("source or target vertex not found: source " + source + ": " + sourceVertex + ", target " + target + ": " + targetVertex);
		}

		insertLink(sourceVertex, targetVertex, type, source);
	}

	private void insertLink(GraphVertex sourceVertex, GraphVertex targetVertex, Link type, DataBean sourceDataBean) {
		if (type.equals(DataBean.Link.GROUPING)) {

			if (sourceVertex.getGroup() != null && targetVertex.getGroup() != null) {
				// both are grouped
				if (sourceVertex.getGroup() == targetVertex.getGroup()) {
					// already in the same group, ignore
					return;
				} else {
					throw new IllegalArgumentException("beans already have different groups");
				}
			}

			if (sourceVertex.getGroup() == null && targetVertex.getGroup() == null) {
				createGroup(sourceVertex.getData()); // create group for the source and add target next
			}

			// other is grouped, other not => group the ungrouped one
			if (sourceVertex.getGroup() != null) {
				sourceVertex.getGroup().addChildVertex(targetVertex);
			} else {
				targetVertex.getGroup().addChildVertex(sourceVertex);
			}

		} else if (type.equals(DataBean.Link.ANNOTATION) || type.equals(DataBean.Link.DERIVATION) || type.equals(DataBean.Link.MODIFICATION)) {

			DefaultEdge linkEdge = new NoLabelEdge(type);
			linkEdge.setSource(sourceVertex.getChildAt(0));
			linkEdge.setTarget(targetVertex.getChildAt(0));

			switch (type) {
			case ANNOTATION:
				setAnnotationEdgeStyle(linkEdge);
				break;
			case DERIVATION:
				setDerivationEdgeStyle(linkEdge);
				break;
			case MODIFICATION:
				setModificationEdgeStyle(linkEdge);
				break;
			}

			this.getGraphLayoutCache().insert(linkEdge);

			if (type.equals(DataBean.Link.DERIVATION) || type.equals(DataBean.Link.MODIFICATION)) {
				layoutManager.updateLayout(sourceVertex, sourceDataBean); // update position if this was made child of other bean
				graphPanel.autoZoom();
				scrollCellToVisibleAnimated(sourceVertex);
				repaint();

			} else if (type.equals(DataBean.Link.ANNOTATION)) {
				moveCloseToAnnotated(sourceVertex.getData());
			}

		} else {
			throw new IllegalArgumentException("unsupported link type: " + type);
		}
	}

	public void removeLink(DataBean source, DataBean target, Link type) {
		GraphVertex sourceVertex = vertexMap.get(source);
		GraphVertex targetVertex = vertexMap.get(target);

		if (type.equals(Link.GROUPING)) {
			sourceVertex.getGroup().removeChildVertex(sourceVertex);

		} else {
			for (DefaultEdge edge : getAllEdgesOfVertex(sourceVertex, this)) {

				// Get link type, source vertex and target vertex
				Link edgeType = (Link) edge.getUserObject();
				GraphVertex edgeSource = (GraphVertex) ((DefaultPort) edge.getSource()).getParent();
				GraphVertex edgeTarget = (GraphVertex) ((DefaultPort) edge.getTarget()).getParent();

				logger.debug("Edge type: " + edgeType + ", edgeSource: " + edgeSource + ", edgeTarget: " + edgeTarget);

				if (edgeSource.equals(sourceVertex) && edgeTarget.equals(targetVertex) && edgeType.equals(type)) {
					// Remove the edge if target, source and link type matched
					graphLayoutCache.remove(new Object[] { edge });
				}
			}
		}
	}

	/**
	 * Gets all edges (incoming and outgoing) of the given vertex
	 * 
	 * @param sourceVertex
	 *            vertex
	 * @param graph
	 *            graph
	 * @return all incoming and outgoing edges of the given vertex
	 */
	private static List<DefaultEdge> getAllEdgesOfVertex(GraphVertex sourceVertex, JGraph graph) {
		
		GraphLayoutCache layoutCache = graph.getGraphLayoutCache();

		List<DefaultEdge> allEdges = new ArrayList<DefaultEdge>();

		// Get outgoing vertexes
		for (Object edge : layoutCache.getOutgoingEdges(sourceVertex, null, false, false)) {
			if (edge instanceof DefaultEdge) {
				allEdges.add((DefaultEdge) edge);
			}
		}

		// Get incoming vertexes
		for (Object edge : layoutCache.getIncomingEdges(sourceVertex, null, false, false)) {
			if (edge instanceof DefaultEdge) {
				allEdges.add((DefaultEdge) edge);
			}
		}

		// Return the list
		return allEdges;
	}

	/**
	 * A normal edge without label. I didn't find any better way to disable the label of the edges...
	 */
	private class NoLabelEdge extends DefaultEdge {

		public NoLabelEdge(Object userObject) {
			super(userObject);
		}

		@Override
		public String toString() {
			return "";
		}
	}

	/*
	 * JGraph isn't able to handle L&F changes when it has content, at least for now. Overriding this method we simply prevent any L&F
	 * changes and the initial L&F will be set in the contructor with a direct call to parents method. This isn't very nice solution, but
	 * hopefully it can be removed later with newer JGraph versions.
	 * 
	 * @see org.jgraph.JGraph#updateUI()
	 */
	@Override
	public void updateUI() {
	}

	@Override
	public void graphChanged(GraphModelEvent e) {
		// listen for vertex moves to store the new position in DataBean
		for (Object o : e.getChange().getChanged()) {
			if (o instanceof GraphVertex) {
				GraphVertex vertex = (GraphVertex) o;
				DataBean bean = getBean(vertex);
				if (bean != null) {					
					bean.setPosition(vertex.getX(), vertex.getY());
				}
			}
		}
	}

	private DataBean getBean(GraphVertex vertex) {
		for (DataBean bean : vertexMap.keySet()) {
			if (vertexMap.get(bean) == vertex) {
				return bean;
			}
		}
		return null;
	}
}
