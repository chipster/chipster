package fi.csc.microarray.client.dataviews.vertexes;

import java.util.ArrayList;
import java.util.List;

import org.jgraph.graph.GraphConstants;

import fi.csc.microarray.client.dataview.MicroarrayGraph;
import fi.csc.microarray.databeans.DataBean;

/**
 * Graph vertex for phenodata. Basicly the PhenodataVertex is just like 
 * GraphVertex, but this class has same extra method compared to GraphVertex 
 * class, and some methods have been taken out of use. The methods which are 
 * not in use (addChildVertex, removeChildVertex) throw an 
 * UnsupportOperationException if they are tried to use.
 * 
 * Maybe the most important method of the class is <code>isPhenodataSet</code>
 * which is used to define wheater the warning sign icon is drawn or not.
 * 
 * Phenodata has also an oval shape and yellow warning sign which are drawn 
 * in the renderer class.
 * 
 * @see GraphRenderer
 * @author mkoski
 *
 */
public class PhenodataVertex extends GraphVertex {

	public PhenodataVertex(int x, int y, DataBean data, MicroarrayGraph graph) {
		super(x, y, data, graph);
		GraphConstants.setSelectable(this.getAttributes(), true);
	}
	
	/** 
	 * Returns "phe".
	 * 
	 * @return "phe"
	 */
	@Override
	public String toString() {
		return "phe";
	}
	
	/**
	 * This method can be called, but it will not do anything.
	 * Phenodata verteces don't have children. 
	 */
	public void addChildVertex(GraphVertex child){
		// empty method
	}

	/**
	 * This method can be called, but it will not do anything.
	 * Phenodata verteces don't have children.
	 * 
	 * @return always false 
	 */
	public boolean removeChildVertex(GraphVertex child){
		return false;
	}

	/**
	 * Returns an empty list. Phenodata verteces don't 
	 * have children. 
	 * 
	 * @return empty list
	 */
	public List<GraphVertex> getChildVertices() {
		return new ArrayList<GraphVertex>();
	}

	/**
	 * Returns an empty list. Phenodata verteces don't 
	 * have children. 
	 * 
	 * @return empty list
	 */
	public List<GraphVertex> getRootVertixesOfThisVertex(){
		return new ArrayList<GraphVertex>();
	}

	/**
	 * Returns whether the phenodata is set or not. The warning sign will became 
	 * visible depending on the return value.
	 */
	public boolean isPhenodataSet(){
		return getData().queryFeatures("/phenodata/is-complete").exists();
	}

	@Override
	public String getToolTipString() {
		if(isPhenodataSet()){
			return super.getToolTipString();
		} else {
			return "Phenodata needs to be completed";
		}
		
	}
}
