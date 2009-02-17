package fi.csc.microarray.analyser.ws.resultgraph;

import java.util.LinkedList;

public class ResultGraphNode {
	
	private LinkedList<ResultGraphNode> children = new LinkedList<ResultGraphNode>();
	private ResultGraphNode parent = null;
	private String type;
	private String name;
	private String id;
	
	public ResultGraphNode(String type, String name, String id) {
		this.type = type;
		this.name = name;
		this.id = id;
	}

	public String getType() {
		return type;
	}

	public String getName() {
		return name;
	}

	public String getId() {
		return id;
	}

	public void addChild(ResultGraphNode child) {
		if (child.parent != null) {
			throw new IllegalArgumentException("child already has parent");
		}
		child.parent = this;
		this.children.add(child);
	}
	
	public void addChildren(Iterable<ResultGraphNode> children) {
		for (ResultGraphNode child : children) {
			addChild(child);
		}
	}
	

}
