package fi.csc.microarray.client.operation;

import java.awt.GridLayout;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import fi.csc.microarray.constants.VisualConstants;

/**
 * A panel that lists operations filtered by a certain criterion.
 * 
 * @author naktinis
 *
 */
@SuppressWarnings("serial")
public class OperationFilterPanel extends JPanel
                                  implements ListSelectionListener {
    
    private static final OperationCategory CATEGORY_ALL = new OperationCategory("All");
    private OperationCategory CATEGORY_NO_MATCHES = new OperationCategory("No matches");
    
    private OperationPanel operationPanel;
    private List<OperationCategory> categories;
    
    private List<OperationDefinition> matchingOperations = new LinkedList<OperationDefinition>();
    private JList categoryList;
    private JList visibleOperationsList;
    
    private ExecutionItem selectedOperation;
    
    public OperationFilterPanel(OperationPanel parent,
           List<OperationCategory> categories) {
        super(new GridLayout(1, 2));
        this.operationPanel = parent;
        this.categories = categories;
        
        categoryList = new JList();
        categoryList.setSelectedIndex(0);
        categoryList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        categoryList.addListSelectionListener(this);
        categoryList.setCellRenderer(new OperationChoicePanel.CategoryListRenderer());
        categoryList.getInsets().right = 1;
        categoryList.setName("categoryList");
        
        // Operation list shows the results
        visibleOperationsList = new JList();
        visibleOperationsList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        visibleOperationsList.addListSelectionListener(this);
        visibleOperationsList.setCellRenderer(new OperationChoicePanel.FontSizeFriendlyListRenderer());
        visibleOperationsList.getInsets().right = 1;
        visibleOperationsList.setName("operationList");
        
        JScrollPane categoryListScroller = new JScrollPane(categoryList);       
        JScrollPane operationListScroller = new JScrollPane(visibleOperationsList);
        categoryListScroller.setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1,
                VisualConstants.OPERATION_LIST_BORDER_COLOR));
        operationListScroller.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));

        this.add(categoryListScroller);
        this.add(operationListScroller);        
    }
    
    /**
     * Load filtered operations.
     */
    public void loadFilteredOperations(String filterPhrase) {
        Vector<OperationDefinition> filteredOperations =
            new Vector<OperationDefinition>();
        LinkedList<Float> filteredOperationsWeights = new LinkedList<Float>();
        
        for (OperationCategory category : categories) {
            selectedOperation = null;
            for (OperationDefinition operation : category.getOperationList()) {
                int indexInTitle = operation.getDisplayName().toLowerCase().
                                   indexOf(filterPhrase.toLowerCase());
                int indexInDescription = operation.getDescription().toLowerCase().
                                   indexOf(filterPhrase.toLowerCase());
                // Phrase found in title
                float weightTitle = (float) Math.min(indexInTitle + 1, 1) /
                    // Favour the beginning of the string
                    (float) (indexInTitle + 2) /
                    // Favour shorter matches
                    (float) operation.getDisplayName().length() * 2;
                // Phrase found in description
                float weightDescription = (float) Math.min(indexInDescription + 1, 1) /
                    (float) operation.getDescription().length();
                
                // Filter operations with positive weights
                if (weightTitle + weightDescription > 0) {
                    filteredOperations.add(operation);
                    filteredOperationsWeights.add(weightTitle + weightDescription);
                }
            }          
        }
        
        // List comparator
        class OperationFilterComparator implements
              Comparator<OperationDefinition> {
            Vector<OperationDefinition> operations;
            LinkedList<Float> weights;
            public OperationFilterComparator(
                    Vector<OperationDefinition> filteredOperations,
                    LinkedList<Float> filteredOperationsWeights) {
                this.operations = filteredOperations;
                this.weights = filteredOperationsWeights;
            }
            
            public int compare(OperationDefinition o1,
                               OperationDefinition o2) {
                float weight1 = weights.get(operations.indexOf(o1));
                float weight2 = weights.get(operations.indexOf(o2));
                return Math.round(Math.signum(weight2 - weight1));
            }
        }
        
        // Sort according to weights
        Collections.sort(filteredOperations,
                new OperationFilterComparator(filteredOperations,
                filteredOperationsWeights));
        
        // Update category list
        LinkedList<OperationCategory> matchingCategories = new LinkedList<OperationCategory>();
        for (OperationDefinition filteredDefinition : filteredOperations) {
        	OperationCategory oc = filteredDefinition.getCategory();
        	if (!matchingCategories.contains(oc)) {
        		matchingCategories.add(oc);
        	}
        }

        // sort categories
        Collections.sort(matchingCategories, new Comparator<OperationCategory>() {

			@Override
			public int compare(OperationCategory category1, OperationCategory category2) {
				return categories.indexOf(category1) - categories.indexOf(category2);
			}
        	
        });

        // add extra categories to the top
        if (filteredOperations.size() > 0) {
        	matchingCategories.addFirst(CATEGORY_ALL);
        } else {
        	CATEGORY_NO_MATCHES.setName("No matching tools for \'" + filterPhrase +"\'");
        	matchingCategories.addFirst(CATEGORY_NO_MATCHES);
        }

        // update visible categories
        categoryList.setListData(matchingCategories.toArray());
        
        // update visible tools
        this.matchingOperations = filteredOperations;
        visibleOperationsList.setListData(filteredOperations);
        
        // select best match
        if (!filteredOperations.isEmpty()) {
        	visibleOperationsList.setSelectedValue(filteredOperations.firstElement(), true);
        }
    }


    /**
     * User has selected an operation.
     */
    public void valueChanged(ListSelectionEvent event) {
        
    	// category selected
    	if (event.getSource() == categoryList && categoryList.getSelectedValue() != null) {
        	OperationCategory selectedCategory = (OperationCategory) categoryList.getSelectedValue();

        	// get matching tools in the selected category
    		List<OperationDefinition> matchingToolsInSelectedCategory = new LinkedList<OperationDefinition>();

    		if (selectedCategory == CATEGORY_NO_MATCHES) {
    			// don't allow selecting this
    			categoryList.clearSelection();
    			return;
    		}
    		
        	// all category
    		else if (selectedCategory == CATEGORY_ALL) {
        		matchingToolsInSelectedCategory = matchingOperations;
        	} 
        	
        	// normal category
        	else {
        		for (OperationDefinition tool : matchingOperations) {
        			if (tool.getCategory().equals(selectedCategory)) {
        				matchingToolsInSelectedCategory.add(tool);
        			}
        		}
        	}

        	// make the matching tools visible
        	visibleOperationsList.setListData(matchingToolsInSelectedCategory.toArray());
        	
        	// select the best match
        	if (!matchingToolsInSelectedCategory.isEmpty()) {
        		visibleOperationsList.setSelectedValue(matchingToolsInSelectedCategory.get(0), true);
        	}

        // tool selected	
        } else if (event.getSource() == visibleOperationsList) {
        	Object selected = visibleOperationsList.getSelectedValue();
        	if (selected instanceof ExecutionItem) {
        		
        		selectedOperation = (ExecutionItem) selected;
        		operationPanel.selectOperation(selectedOperation);
        		
        		// also select the category of this operation
        		// avoid category changed event which would cause limiting the visible
        		// tools to this category only
        		categoryList.removeListSelectionListener(this);
        		categoryList.setSelectedValue(((OperationDefinition)selected).getCategory(), true);
        		categoryList.addListSelectionListener(this);
        	}
        }
    }
}
