package fi.csc.microarray.client.dataview;

import java.awt.Cursor;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;

import javax.swing.SwingUtilities;

import org.apache.log4j.Logger;
import org.jgraph.JGraph;
import org.jgraph.graph.BasicMarqueeHandler;
import org.jgraph.graph.CellView;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataviews.vertexes.AbstractGraphVertex;
import fi.csc.microarray.client.dataviews.vertexes.GraphVertex;
import fi.csc.microarray.client.dataviews.vertexes.GroupVertex;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.databeans.DataBean;

/**
 * This class takes care of the selection made on the graph vertex. Basicly there 
 * are two types of selections: normal selections (mouse click over the cell) and 
 * marquee selections (mouse pressed and dragged and all the cells inside the 
 * rectangle are selected). 
 * 
 * There is few things in this class which may not be very simple. The important method 
 * of this class is <code>isForceMarqueeTigger</code> which gains mouse event control 
 * to this class if it return <code>true</code>. This method prevents complications 
 * between this and JGraph's own mouse handlers.
 * 
 * Other important method is <code>isMarqueeTiggerEvent</code> which defines whether 
 * the mouse press event is a marquee selection (the rectangle is drawn) or not.
 * 
 * When the marquee selection is made the <code>handleMarqueeEvent</code> method is 
 * called.
 * 
 * After the selection is made the Application class fires an DatasetChoice event 
 * which is listened and graph's cells' borders are drawn thiner or thicker if 
 * they are selected or not.
 * 
 * If user clicks on an empty space the click is registered as an empty marquee 
 * event and it is handeled in the <code>handleMarqueeEvent</code> method
 * 
 * To attach this tool it must be registered graph's MarqueeSelectionListener 
 * and MouseListener. This is done by graph's addMarqueeListener and 
 * addMouseListener methods.
 * 
 * @author klemela, mkoski
 *
 */
public class SelectionTool extends BasicMarqueeHandler implements MouseListener{
	
	private MicroarrayGraph graph;
	
	private static final Logger logger = Logger.getLogger(SelectionTool.class);
	private ClientApplication application = Session.getSession().getApplication();
	
	/**
	 * Constructor
	 * 
	 * @param graph Used graph
	 */
	public SelectionTool(MicroarrayGraph graph) {
		this.graph = graph;
	}
	
	/**
	 * Overrides parent method and casts JGraph to MicroarrayGraph
	 */
	public void handleMarqueeEvent(MouseEvent e, JGraph graph,
			Rectangle2D bounds) {
		this.handleMarqueeEvent(e, (MicroarrayGraph)graph, bounds);
	}
	
	/**
	 * This method is called when marquee selection is made (mouse is released). 
	 * It gets the list of all vertexes inside the selection rectangle and 
	 * passes it the to <code>setSelected</code> or <code>addSelected</code> method
	 * depending on control key state.
	 * 
	 * @param e 		Mouse event
	 * @param graph 	Graph
	 * @param bounds 	Selection rectangle
	 */
	public void handleMarqueeEvent(MouseEvent e, MicroarrayGraph graph,
			Rectangle2D bounds) {
				
		logger.debug("handleMarqueeEvent: isControlDown " + e.isControlDown());
		
		List<AbstractGraphVertex> vertexesToSelect = new ArrayList<AbstractGraphVertex>(); 
		
		// Iterates all cell views (cells that are visible)
		CellView[] views = graph.getGraphLayoutCache().getCellViews();
		for (CellView view : views) {
			if(view != null){
				if(view.getCell() instanceof GroupVertex){
					GroupVertex group = (GroupVertex)view.getCell();
					
					// Checks that group vertex is in marquee bounds and is collapsed
					// If group is not collapsed the selection is made like normal vertexes
					if( bounds.contains(group.getBounds())){
						vertexesToSelect.add(group);
						// this.selectGroup(group, e);
					} else {
						// Ignore
					}
				} 
				
				else if(view.getCell() instanceof AbstractGraphVertex){
					AbstractGraphVertex vertex = (AbstractGraphVertex)view.getCell();
					
					// Checks that vertex is in marquee bounds
					if(bounds.contains(vertex.getBounds())){
						vertexesToSelect.add(vertex);
					}
				}
			}
		}
		
		// Passes the list of vertexes inside the selection rectangle to select method
		if(e.isControlDown()){
			this.setSelected(vertexesToSelect);
		} else {
			this.addSelected(vertexesToSelect);
		}
	}
	
	/**
	 * Selects a single cell. Actually this method only creates a list and adds 
	 * the cell to the it and passes the list to the other select method.
	 * 
	 * @param vertex 	cell to select
	 * @param e 		mouse event
	 */
	private void setSelected(AbstractGraphVertex vertex){
		List<AbstractGraphVertex> vertexesToSelect = new ArrayList<AbstractGraphVertex>();
		vertexesToSelect.add(vertex);
		this.setSelected(vertexesToSelect);
	}
	
	private boolean isSelected(AbstractGraphVertex vertex){
		if(vertex instanceof GraphVertex){
			return application.getSelectionManager().getSelectedDataBeans().contains(((GraphVertex)vertex).getData());
		} else {
			return false;
		}
	}
	
	/**
	 * Sets the cells as selected cells. This method clears all the previous 
	 * selections user has made.
	 */
	private void setSelected(List<AbstractGraphVertex> vertexesToSelect){
		
		// Sets selection mode
		boolean multiSelectionMode;
		if (vertexesToSelect.size() > 1 || (vertexesToSelect.size() == 1 &&	vertexesToSelect.get(0) instanceof GroupVertex)) {
			// User selected more than one vertex or a group
			multiSelectionMode = true;
		} else {
			// User selected only one vertex and is not pressing control key
			multiSelectionMode = false;
		}
		
		// First clear all previous selections
		application.getSelectionManager().clearAll(false, graph);
		
		// Then do the selection
		select(vertexesToSelect, multiSelectionMode);
	}
	
	/**
	 * Sets the passed group as selected. This method clears all the previous 
	 * selections user has made.
	 */
	private void setSelectedGroup(GroupVertex group, MouseEvent e){
		setSelected(group.getChildVertexes());
	}
	
	/**
	 * Adds vertexes to selected. Does not clear previous selections
	 * @param vertexesToSelect
	 * @param e
	 */
	private void addSelected(List<AbstractGraphVertex> vertexesToSelect){
		select(vertexesToSelect, true);
	}
	
	/**
	 * Adds vertex to selected. Does not clear previous selections
	 * 
	 * @param vertexesToSelect
	 * @param e
	 */
	private void addSelected(AbstractGraphVertex vertex){
		List<AbstractGraphVertex> vertexesToSelect = new ArrayList<AbstractGraphVertex>();
		vertexesToSelect.add(vertex);
		this.addSelected(vertexesToSelect);
	}
	
	/**
	 * Selects the given vertexes using the given selection mode. This 
	 * method shouldn't be used directly. The setSelected and addSelected 
	 * methods should be use this method.
	 * 
	 * @param vertexesToSelect
	 * @param multiSelectionMode
	 */
	private void select(List<AbstractGraphVertex> vertexesToSelect, boolean multiMode) {
		if (!multiMode) {
			application.getSelectionManager().clearAll(false, graph);
		}
		
		for (AbstractGraphVertex vertex : vertexesToSelect) {
			
			if (vertex instanceof GroupVertex) {
				for (DataBean selectedItem : ((GroupVertex)vertex).getChildrenData()) {
					application.getSelectionManager().selectMultiple(selectedItem, graph);
				}
				
			} else {
				DataBean data = ((GraphVertex)vertex).getData();
				if(application.getSelectionManager().isSelected(data)){
					application.getSelectionManager().deselectMultiple(data, graph);
				} else {
					application.getSelectionManager().selectMultiple(data, graph);
				}
			}
		}
	}
	
	/**
	 * Gains full control of mouse event prosessing. Every mouse event is forced 
	 * to be an marquee event which means that all the mouse event processing 
	 * is done in this class.
	 * 
	 * If the left mouse is clicked and it is not a marquee selection and is not 
	 * in the folding handle hit region the method 
	 * will return <code>false</code>, which means that JGraph BasicUI handles the 
	 * event in its own way. This way the cells are moveable by mouse dragging.
	 * 
	 */
	@Override
	public boolean isForceMarqueeEvent(MouseEvent e) {
		if(SwingUtilities.isLeftMouseButton(e) && !isMarqueeTriggerEvent(e, graph)){
			return false;
		} else {
			return true;
		}
	}
	
	/**
	 * Overrides isMarqueeTiggerEvent method so that it return false if 
	 * user is trying to start marquee selection over a dataset or user 
	 * clicked group's collapse handle. In this case 
	 * marquee selection is disabled and dataset vertex is moveable or group is collapsed.
	 */
	@Override
	public boolean isMarqueeTriggerEvent(MouseEvent e, JGraph graph) {
		
		Object clickedCell = graph.getFirstCellForLocation(e.getX(), e.getY());
		
		// User clicked a group or dataset
		if(clickedCell != null){
			return false;
		} 
		
		// User clicked empty space. It's ok to start marquee
		else {
			return super.isMarqueeTriggerEvent(e, graph);
		}
	}
	
	/**
	 * Keeps the mouse cursor as a default cursor if right mouse button is clicked or 
	 * user has made a double click. 
	 * 
	 */
	@Override
	public void mousePressed(MouseEvent e) {
		
		// Left mouse button pressed and it is an marquee
		if(SwingUtilities.isLeftMouseButton(e) && isMarqueeTriggerEvent(e, graph)){
			// Before marquee rectangle is drawn clear the previous selections. 
			// This also clears the selections if user just clicks to the empty space 
			// on the graph
			application.getSelectionManager().clearAll(true, graph);
			super.mousePressed(e);
		}
		
		// Not a marquee
		else if(SwingUtilities.isRightMouseButton(e)){
			// Change the cursor back from marquee cross
			this.graph.setCursor(Cursor.getDefaultCursor());
			
			// Select the graph under the cursor if there is no datasets 
			// selected already
			this.rightMouseButtonPressed(e);
		}
	}
	
	/**
	 * Consumes the event if is not a marquee event. This prevent the 
	 * BasicGraphUI.MouseHandler to catch the mouse event.
	 */
	@Override
	public void mouseReleased(MouseEvent e) {
		if(SwingUtilities.isLeftMouseButton(e) && isMarqueeTriggerEvent(e, graph)){
			super.mouseReleased(e);
		}

		else {
			e.consume();
		}
		this.graph.setCursor(Cursor.getDefaultCursor());
	}
	
	/**
	 * Prevents the selection rectangle to appear if mouse is dragged 
	 * pressing right mouse button
	 */
	@Override
	public void mouseDragged(MouseEvent e) {
		if(SwingUtilities.isLeftMouseButton(e)){
			super.mouseDragged(e);
		}
	}
	
	/**
	 * Mouse button is clicked, selects the dataset if it is not already selected. 
	 * 
	 * There is some quite complicated checkings that must be done for the selection. The 
	 * main problem is that expanded groups should be not selected at all but for some 
	 * reason <code>getFirstCellForLocation</code> returns expanded groups too. That's 
	 * why the getFirst and getNextCellForLocation used.
	 */
	private void mouseButtonClicked(MouseEvent e){
		// Get cell
		Object cell = graph.getFirstCellForLocation(e.getX(), e.getY());
		
		AbstractGraphVertex vertex = null;
		if (cell instanceof AbstractGraphVertex){
			vertex = (AbstractGraphVertex)cell;
		}
			
		// Okey! Let's select the vertex if it's not selected already
		if(vertex != null){
			if(e.isControlDown()){
				addSelected(vertex);
			} else {
				setSelected(vertex);
			}
		}
	}
	
	private void rightMouseButtonPressed(MouseEvent e){
		// Get cell
		Object cell = graph.getFirstCellForLocation(e.getX(), e.getY());
		
		AbstractGraphVertex vertex = null;
		if (cell instanceof AbstractGraphVertex){
			vertex = (AbstractGraphVertex)cell;
		}
			
		// Get data
		if(vertex != null){
			if(e.isControlDown()){
				addSelected(vertex);
			} else {
				if(!isSelected(vertex)){
					setSelected(vertex);
				}
			}
		}
	}
	
	/**
	 * When users makes a double click the clicked dataset is visualized with the 
	 * best method. If user double clicks a group the group is selected but not 
	 * visualized.
	 * 
	 * @param e
	 */
	private void mouseButtonDoubleClicked(MouseEvent e){
		// Get the clicked cell
		Object cell = graph.getFirstCellForLocation(e.getX(), e.getY());
		
		logger.debug("Selected cell: " + cell);
		
		// Do not visualise collapsed group
		if(cell instanceof GraphVertex){
			application.visualiseWithBestMethod(FrameType.MAIN);
		} else if (cell instanceof GroupVertex ) {
			// User double clicked a collapsed group, selects the group but does not visualise
			setSelectedGroup((GroupVertex)cell, e);
		}
	}

	public void mouseClicked(MouseEvent e){	
		
		// Double click
		if(SwingUtilities.isLeftMouseButton(e) && e.getClickCount() > 1){
			// Change the cursor back from marquee cross
			this.graph.setCursor(Cursor.getDefaultCursor());
			this.mouseButtonDoubleClicked(e);
			
//			 if the Left mouse button is clicked is checked after double-click to prevent useless event from the second click 
		} else if(SwingUtilities.isLeftMouseButton(e) && !isMarqueeTriggerEvent(e, graph)){
				mouseButtonClicked(e);
			}
	}

	public void mouseEntered(MouseEvent e) {
		// Do nothing
		
	}

	public void mouseExited(MouseEvent e) {
		// Do nothing
		
	}
}
