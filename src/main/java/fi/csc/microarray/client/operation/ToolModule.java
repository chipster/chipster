package fi.csc.microarray.client.operation;

import java.util.LinkedList;
import java.util.List;

public class ToolModule {

    private List<ToolCategory> visibleCategories = new LinkedList<ToolCategory>();
    public String getModuleName() {
		return moduleName;
	}

	private List<ToolCategory> hiddenCategories = new LinkedList<ToolCategory>();
    
	private String moduleName;

    public ToolModule(String moduleName) {
		this.moduleName = moduleName;
	}

    /**
     * @return true iff module should be shown in GUI 
     */
    public boolean isVisible() {
    	return !visibleCategories.isEmpty();
    }
    
	/**
     * @return categories that are visible to end-user.
     */
    public List<ToolCategory> getVisibleCategories() {
        return visibleCategories;
    }

    /**
     * @return categories that are hidden from end-user, but still
     * available for execution.
     */
    public List<ToolCategory> getHiddenCategories() {
        return hiddenCategories;
    }
    
    public void addVisibleToolCategory(ToolCategory category) {
    	visibleCategories.add(category);
    }

    public void addHiddenToolCategory(ToolCategory category) {
    	hiddenCategories.add(category);
    }

    public OperationDefinition getOperationDefinition(String toolId) {
    	for (List<ToolCategory> categories : new List[] { visibleCategories, hiddenCategories}) {
    		for (ToolCategory category : categories) {
    			for (OperationDefinition tool : category.getToolList()) {
    				if (tool.getID() == toolId) {
    					return tool;
    				}
    			}
    		}
    	}
    	
    	return null;
		
	}

}
