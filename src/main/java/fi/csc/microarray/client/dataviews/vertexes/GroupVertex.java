package fi.csc.microarray.client.dataviews.vertexes;

import java.awt.Point;
import java.util.ArrayList;
import java.util.List;

import org.jgraph.graph.GraphConstants;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataview.MicroarrayGraph;
import fi.csc.microarray.databeans.DataBean;


/**
 * Vertex for a group cell
 * 
 * @author mkoski
 *
 */
public class GroupVertex extends AbstractGraphVertex {
	
	private ClientApplication application = Session.getSession().getApplication();			
	
	public GroupVertex(int x, int y, GraphVertex[] children, MicroarrayGraph graph) {
		super(x, y, null, graph);
		
		this.addPort();
		
		// Set group vertex to back
		graph.getGraphLayoutCache().toBack(new Object[] { this });
		
		if(children != null){
			for(GraphVertex vertex : children){
				addChildVertex(vertex);
			}
		}
		
		graph.repaint();
	}
	
		
	public void collapse(){
		this.getGraph().getGraphLayoutCache().collapse(new Object[] { this });		
		getGraph().updateSelectedCells();
	}
	
	/* There isn't anymore need to show expanded group, but it has to be opened
	 * when new vertexes are added, because adding vertexes to closed group
	 * causes strange things.
	 */	
	public void expand(){
		this.getGraph().getGraphLayoutCache().expand(new Object[] { this });
		getGraph().updateSelectedCells();
	}
	
	/**
	 * Gets the child vertexes of this group.
	 * 
	 * @return child vertexes
	 */
	public List<AbstractGraphVertex> getChildVertexes(){
		List<AbstractGraphVertex> vertexes = new ArrayList<AbstractGraphVertex>();
		List children = this.getChildren();
		for (Object child : children) {
			if(child instanceof AbstractGraphVertex){
				vertexes.add((AbstractGraphVertex)child);
			}
		}
		return vertexes;
		
	}

	public List<DataBean> getChildrenData() {
		List<DataBean> childData = new ArrayList<DataBean>();
		for (Object child : this.getChildren()) {
			if(child instanceof GraphVertex){
				GraphVertex vertex = (GraphVertex)child;
				assert(vertex.getData() instanceof DataBean);
				childData.add((DataBean)vertex.getData());
			}
		}
		return childData;
	}
		
	@Override
	public String toString(){
		
		// For some reason child count is too big by one number
		int childCount = getChildCount() - 1;
		
		if(childCount > 1){
			return childCount + " files";
		} else {
			return childCount + " file";
		}
	}

	@Override
	public String getToolTipString() {
		return this.toString();
	}
	
	public boolean isSelectable(){
		return (Boolean)this.getAttributes().get(GraphConstants.SELECTABLE);
	}
	
	/**
	 * Adds vertex to group. There is no list of group members. The group 
	 * members are stored to jgraph's parent maps.
	 */
	public void addChildVertex(GraphVertex vertex) {
		
		// Expand the group before adding the vertex, otherwise the
		// group vertex appears in strange place
		this.expand();		
		
		//Group seems to move with it's children, so we give the same location for them 
		vertex.setPosition(new Point((int)this.getBounds().getX(), (int)this.getBounds().getY()));

		this.add(vertex);
				
		this.collapse();		
	}
	
	@Override
	public void setPosition(Point p){
		for(AbstractGraphVertex child : getChildVertexes()){
			child.setPosition(p);
		}
		super.setPosition(p);
	}
		
	public boolean isAllChildrenSelected(){
		List<DataBean> selectedDatas = application.getSelectionManager().getSelectedDataBeans();
		
		// No selected data at all
		if(selectedDatas.size() == 0){
			return false;
		}
		
		for (DataBean childData : this.getChildrenData()) {
			if(!selectedDatas.contains(childData)){
				return false;
			}
		}
		return true;
	}

	@Override
	public boolean removeChildVertex(GraphVertex vertex) {
		if(!this.contains(vertex)){
			return false;
		}
		
		this.expand();
		
		if(getChildVertexes().size() == 0){
			this.getGraph().getGraphLayoutCache().remove(new Object[] { this } );
		}
		this.collapse();
		return true;
	}
	
	public boolean contains(GraphVertex vertex){
		return vertex.getGroup().equals(this);
	}


	@Override
	public int getDefaultHeight() {
		return DEFAULT_HEIGHT; 
	}


	@Override
	public int getDefaultWidth() {
		return DEFAULT_WIDTH;	
	}
}
