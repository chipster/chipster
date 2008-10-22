package fi.csc.microarray.proto.repository.schema;

import java.util.ArrayList;
import java.util.List;

public abstract class ParameterItem {
	
	private String name;
	private String description;
	private ParameterClass parent;
	
	public ParameterItem(String name) {
		this.name = name;
		this.parent = null;
	}
	
	public String getName() {
		return name;
	}
	
	public String getDescription() {
		return description;
	}
	
	public void setDescription(String description) {
		this.description = description;
	}
	
	public void setParent(ParameterClass parent) {
		this.parent = parent;
	}
	
	public ParameterClass getParent() {
		return parent;
	}
	
	public String toString() {
		return name;
	}
	
	public List<ParameterItem> getChoiceList() {
		List<ParameterItem> list;
		if (this.parent != null) {
			list = this.parent.getChoiceList();
		} else {
			list = new ArrayList<ParameterItem>();
		}
		list.add(this);
		return list;
	}
	
	public List<List<ParameterItem>> getPossiblePathsToNamedDescendant(String entityName) {
		if (this.getName().toLowerCase().indexOf(entityName.toLowerCase()) != -1) {
			// If this is what we're looking for, create the new list:
			List<ParameterItem> pathList =
				new ArrayList<ParameterItem>();
			pathList.add(this);
            // Then create the list of lists containing the list we just created:
            List<List<ParameterItem>> listOfPathLists =
                new ArrayList<List<ParameterItem>>();
            listOfPathLists.add(pathList);
			return listOfPathLists;
		} else {
			if (this instanceof ParameterClass) {
                List<List<ParameterItem>> currentListOfPathLists = null;
				if (((ParameterClass) this).getInstances() != null) {
					for (ParameterInstance inst : ((ParameterClass) this).getInstances()) {
						List<List<ParameterItem>> listOfPathLists =
							inst.getPossiblePathsToNamedDescendant(entityName);
						if (listOfPathLists != null) {
                            if (currentListOfPathLists == null) {
                                currentListOfPathLists = listOfPathLists;
                            } else {
                                currentListOfPathLists.addAll(listOfPathLists);
                            }
						}
					}
				}
				if (((ParameterClass) this).getSubclasses() != null) {
                    for (ParameterClass c : ((ParameterClass) this).getSubclasses()) {
                        List<List<ParameterItem>> listOfPathLists =
                            c.getPossiblePathsToNamedDescendant(entityName);
                        if (listOfPathLists != null) {
                            if (currentListOfPathLists == null) {
                                currentListOfPathLists = listOfPathLists;
                            } else {
                                currentListOfPathLists.addAll(listOfPathLists);
                            }
                        }
                    }  
				}
                
                if (currentListOfPathLists != null) {
                    for (List<ParameterItem> list : currentListOfPathLists) {
                        list.add(0, this);
                    }
                }
                
				// Return the search results, augmented with this node:
                return currentListOfPathLists;
			}
			return null;
		}
	}
}