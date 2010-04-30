package fi.csc.microarray.client.operation;

import java.awt.GridLayout;
import java.util.Collection;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.apache.log4j.Logger;

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
        
        for (OperationCategory category : categories) {
            selectedOperation = null;
            for (OperationDefinition operation : category.getOperationList()) {
                if (operation.getID().indexOf(filterPhrase) != -1) {
                    filteredOperations.add(operation);
                }
            }
            
            parent.enableAction(false);
        }
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
