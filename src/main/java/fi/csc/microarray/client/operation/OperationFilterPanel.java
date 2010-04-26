package fi.csc.microarray.client.operation;

import java.awt.GridLayout;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.apache.log4j.Logger;

import sun.awt.image.IntegerComponentRaster;

/**
 * A panel that lists operations filtered by a certain criterion.
 * 
 * @author naktinis
 *
 */
@SuppressWarnings("serial")
public class OperationFilterPanel extends JPanel
                                  implements ListSelectionListener {
    
    // Logger for this class
    private static final Logger logger = Logger
        .getLogger(OperationChoicePanel.class);
    
    private OperationPanel parent;
    private Collection<OperationCategory> categories;
    
    private JList operationList;
    
    private ExecutionItem selectedOperation;
    
    public OperationFilterPanel(OperationPanel parent,
           Collection<OperationCategory> categories) {
        super(new GridLayout(1, 1));
        this.parent = parent;
        this.categories = categories;
        
        operationList = new JList();
        operationList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
        operationList.addListSelectionListener(this);
        operationList.setCellRenderer(new OperationChoicePanel.FontSizeFriendlyListRenderer());
        //operationList.addMouseListener(new MouseClickListener());
        operationList.getInsets().right = 1;
        operationList.setName("operationList");
        
        JScrollPane operationListScroller = new JScrollPane(operationList);
        operationListScroller.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
        
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
                int indexInTitle = operation.getName().indexOf(filterPhrase);
                int indexInDescription = operation.getDescription().indexOf(filterPhrase);
                // Phrase found in title
                float weightTitle = (float) Math.min(indexInTitle + 1, 1) /
                    // Favour the beginning of the string
                    (float) (indexInTitle + 2) /
                    // Favour shorter matches
                    (float) operation.getName().length() * 2;
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
        
        parent.enableAction(false);
        
        logger.debug("found " + filteredOperations.size() + " operations " +
                     "for phrase \"" + filterPhrase + "\"");
        operationList.setListData(filteredOperations);
        parent.selectOperation(null);
    }

    /**
     * User has selected an operation.
     */
    public void valueChanged(ListSelectionEvent arg0) {
        Object selected = operationList.getSelectedValue();
        if (selected instanceof ExecutionItem) {
            selectedOperation = (ExecutionItem) selected;
            parent.selectOperation(selectedOperation);
        }
    }
}
