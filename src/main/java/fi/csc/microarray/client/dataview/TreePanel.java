package fi.csc.microarray.client.dataview;

import java.awt.CardLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTree;
import javax.swing.SwingUtilities;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.event.TreeSelectionListener;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.DefaultTreeCellRenderer;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreeNode;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataItemCreatedEvent;
import fi.csc.microarray.databeans.DataItemRemovedEvent;

/**
 * A GUI component for viewing the DataBeans in a tree structure,
 * organized in parent elements called DataFolders. Within each folder,
 * the datasets are kept in alphabetical order.
 * 
 * @author Janne Käki, Aleksi Kallio
 * 
 */
public class TreePanel extends JPanel implements DataChangeListener, TreeSelectionListener, PropertyChangeListener {
	
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(TreePanel.class);
	
    /**
     * A listener class that will notice whenever the user selects a node of
     * the tree. Then it will extract the corresponding DataBean or DataFolderGUIElement
     * and call selectElement method for appropriate action. In the case of
     * DataSets, this method will also notify registered DataSetChoiceListeners.
     */
    private class TreeMouseListener extends MouseAdapter {
        
    	@Override
        public void mousePressed(MouseEvent e) {
    		
    		maybeShowPopup(e);
        	// double click
        	if (SwingUtilities.isLeftMouseButton(e) && e.getClickCount() > 1) {
        		
        		DataItem selectedItem = this.getSelectedElementFrom(e);
        		
        		if (selectedItem instanceof DataBean) {        			
        			application.visualiseWithBestMethod(FrameType.MAIN);
        			
        		} else if (selectedItem instanceof DataFolder) {
        			// select all child beans
        			DataFolder folder = (DataFolder)selectedItem;
        			application.getSelectionManager().clearAll(false, this);
        			application.getSelectionManager().selectMultiple(folder.getChildren(), this);
        		}
        	} 
        }

    	@Override
        public void mouseReleased(MouseEvent e) {
        	maybeShowPopup(e);        
        }
        
        private DataItem getSelectedElementFrom(MouseEvent e) {
            TreePath selectedPath = tree.getPathForLocation(e.getX(), e.getY());
            if (selectedPath != null) {
                Object selectedNode = selectedPath.getLastPathComponent();
                if (selectedNode instanceof DefaultMutableTreeNode) {
                    Object selectedUserObject =
                        ((DefaultMutableTreeNode) selectedNode).getUserObject();
                    if (selectedUserObject != null &&
                            selectedUserObject instanceof DataItem) {
                        return (DataItem) selectedUserObject;
                    }
                }
            }
            return null;
        }
        
        private void maybeShowPopup(MouseEvent e){
        	 if (e.isPopupTrigger()) {
        		 if (getSelectedElementFrom(e) instanceof DataFolder){
        			 //No special handling for the multiple folder selection and we show it
        			 //by removing all other than last selection
        	         tree.setSelectionPath(tree.getPathForLocation(e.getX(), e.getY()));
        			 
         			application.showPopupMenuFor(e, getSelectedElementFrom(e));
         			
         		} else {
         			List<DataItem> items = new ArrayList<DataItem>();
         			for (DataBean bean : application.getSelectionManager().getSelectedDataBeans()) {
         				items.add(bean);
         			}
         			application.showPopupMenuFor(e, items);
         		}
             }
        }
    }
    
    /**
     * Helper class that defines how the tree nodes are to be rendered on
     * screen. Is in fact a JPanel which customizes its appearance (icon,
     * font, and font decorations) according to each dataset's properties.
     */
    private class CustomTreeCellRenderer extends DefaultTreeCellRenderer {
        
        public Component getTreeCellRendererComponent(JTree tree,
                Object value, boolean selected, boolean expanded,
                boolean leaf, int row, boolean hasFocus) {
            
            DefaultMutableTreeNode node = (DefaultMutableTreeNode) value;
            DataItem element = (DataItem) node.getUserObject();

            DefaultTreeCellRenderer renderer = 
            	(DefaultTreeCellRenderer)super.getTreeCellRendererComponent(tree, value, selected, 
            			expanded,leaf, row, hasFocus);
            
            renderer.setIcon(application.getIconFor(element));
            renderer.setText(element.toString() + " ");
           
          return this;
        }
    }
    
	private JTree tree = null;
	private DefaultTreeModel treeModel = null;
	private DefaultMutableTreeNode rootNode = null;
	private DataFolder rootFolder = null;
	private Map<DataItem, DefaultMutableTreeNode> nodeMap = new HashMap<DataItem, DefaultMutableTreeNode>();
	
	private ClientApplication application = Session.getSession().getApplication();	
	
	private boolean disableSelectionReporting = false;
	
	private CardLayout cardLayout;
	private JPanel cardParent;
	
	/**
	 * Creates a new TreePanel component, initially empty with only
	 * the default root folder present.
	 */
	public TreePanel(DataFolder rootFolder, JPanel parent, CardLayout cardLayout) {
		super(new GridBagLayout());
		
		this.cardLayout = cardLayout;
		this.cardParent = parent;
		
		// initialise root and data structures
		this.rootFolder = rootFolder;
		this.tree = getTree();
		
		// Setting name for tests
		this.tree.setName("datasetTree");
				
		// create panels and layout
        this.setBorder(BorderFactory.createEmptyBorder());
        this.setMinimumSize(new Dimension(0,0));
        this.setPreferredSize(new Dimension(VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.TREE_PANEL_HEIGHT));
		JPanel dataTreePanel = new JPanel(new GridBagLayout());
		dataTreePanel.setBackground(tree.getBackground());
		GridBagConstraints c = new GridBagConstraints();
		c.gridy = 0;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.gridy = 1;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 1.0; c.weighty = 1.0;
		dataTreePanel.add(tree, c);
		
		// create scoller
		JScrollPane scroller = new JScrollPane(dataTreePanel);
		// Default scroll speed is very slow
		scroller.getVerticalScrollBar().setUnitIncrement(scroller.getVerticalScrollBar().getUnitIncrement()*4);
		
		scroller.setBorder(BorderFactory.createEmptyBorder());
		scroller.setMinimumSize(new Dimension(0,0));
		this.add(scroller, c);
	
		// start listening
		application.getDataManager().addDataChangeListener(this);
		application.addClientEventListener(this);
	}
	
	public Vector<Component> getFocusComponents(){
		Vector<Component> order = new Vector<Component>();
		order.add(tree);		
		return order;
	}
    
    /**
     * @return The data tree instance variable of this TreePanel
     *         (if previously not present, one is created).
     */
    private JTree getTree() {
        if (tree == null) {
            
            this.rootNode = new DefaultMutableTreeNode(rootFolder);
            this.nodeMap.put(rootFolder, rootNode);
            this.treeModel = new DefaultTreeModel(rootNode);
            this.tree = new JTree(treeModel);

            tree.setCellRenderer(new CustomTreeCellRenderer());
            
            this.tree.getSelectionModel().setSelectionMode(
            		TreeSelectionModel.DISCONTIGUOUS_TREE_SELECTION);
            
            tree.addTreeSelectionListener(this);
            tree.addMouseListener(new TreeMouseListener());                    	
        }
        return tree;
    }
	
	/**
	 * Inserts a new DataItem to this tree.
	 * 
	 * @param data New DataItem to be inserted to this view.
	 * @throws DataInsertionException If the parent of the inserted DataBean
	 * 		   was not already present in this structure. DataSets must be
	 * 		   inserted in such order that parent is added first before any
	 *         of its children.
	 */
	private void insertData(DataItem data) {

		if (nodeMap.containsKey(data)) {
			throw new RuntimeException(data.getName() + " already exists in tree");
		} 
		
		if (data.getParent() == null || nodeMap.containsKey(data.getParent())) {

			if (data instanceof DataBean) {
				DefaultMutableTreeNode node = createNode(data);
				nodeMap.put(data, node);				

			} else if (data instanceof DataFolder) {
				DataFolder folder = (DataFolder) data;
				if (folder.getParent() == null) {
					rootFolder.addChild(folder);
				}
				DefaultMutableTreeNode node = createNode(folder);
				nodeMap.put(folder, node);				
				for (DataItem data1 : folder.getChildren()) {
					insertData(data1);
				}

			} else {
				throw new RuntimeException("illegal data type: " + data.getClass().getSimpleName());
			}

		} else {
			throw new RuntimeException("parent of " + data.getName() + " was not inserted (it is " + data.getParent().getName() + ")");
		}
	}

	/**
	 * Removes the given dataset from the tree.
	 * 
	 * @param data The dataset to be removed.
	 * @return True if deletion succeeded, false if it didn't (because
	 * 		   the corresponding tree node was not found).
	 */
	public void removeData(DataItem data) {
		DefaultMutableTreeNode node = nodeMap.get(data);
		if (node != null) {
			treeModel.removeNodeFromParent(node);
			TreeNode parent = node.getParent();
			treeModel.nodeStructureChanged(parent);
			nodeMap.remove(data);
			tree.repaint();
			
		} else {
			throw new RuntimeException(data.getName() + " does not exists in tree");
		}
	}
	
	/**
	 * Creates a new tree node for the given dataset and places it as
	 * a new child to the already existing node corresponding the parent
	 * dataset of the given data. If parent is null, the node and its
	 * data is added to the root folder of the tree.
	 * 
	 * @param data New dataset to be added to the tree.
	 * @return The newly created node for the given dataset (next to
	 * 		   be added to the nodeMap of this TreePanel).
	 */
	private DefaultMutableTreeNode createNode(DataItem data) {
		DefaultMutableTreeNode node = new DefaultMutableTreeNode(data);
		DataFolder parentFolder = data.getParent();
		
		DefaultMutableTreeNode parentNode;
		if (parentFolder != null) { 			
			parentNode = nodeMap.get(parentFolder); // has parent			
		} else {			
			parentNode = rootNode; // has no parent, use root
		}
		
		// insert into tree model
		assert(parentNode != null);
		treeModel.insertNodeInto(node, parentNode, parentNode.getChildCount());
		
		disableSelectionReporting = true;
		TreePath[] paths = tree.getSelectionPaths();
		
		//Internal implemention removes selection (at least in Java 1.5), which is 
		//prevented with preceding and following two lines
		treeModel.nodeStructureChanged(parentNode);
		
		tree.setSelectionPaths(paths);
		disableSelectionReporting = false;
		
		return node;
	}
	
	public void updateNameFor(DataItem data) {
		DefaultMutableTreeNode node = nodeMap.get(data);
		if (node != null) {
			treeModel.nodeStructureChanged(node);
		}
	}

	public void propertyChange(PropertyChangeEvent dataEvent) {
		
		logger.debug("got " + dataEvent.getClass().getSimpleName());
		
        if (dataEvent instanceof DatasetChoiceEvent && 
        		dataEvent.getSource() != this) {
        	
        	LinkedList<TreePath> paths = new LinkedList<TreePath>();
        	for(DataBean bean : application.getSelectionManager().getSelectedDataBeans()){
        		paths.add(new TreePath(nodeMap.get(bean).getPath()));
        	}        	
        	TreePath path[] = paths.toArray(new TreePath[paths.size()]);
        	
        	disableSelectionReporting  = true;
        	tree.setSelectionPaths(path);
        	disableSelectionReporting = false;
        	
        	if(path.length > 0){
        		tree.scrollPathToVisible(path[path.length - 1]);
        	}
        }        
	}
	
	public void dataChanged(DataChangeEvent dataEvent) {
		
		logger.debug("got " + dataEvent.getClass().getSimpleName());
		
        if (dataEvent instanceof DataChangeEvent) {
        	
        	DataItem data = ((DataChangeEvent)dataEvent).getDataItem();
        	
        	if (dataEvent instanceof DataItemCreatedEvent) {
        		insertData(data);
        		
        	} else if (dataEvent instanceof DataItemRemovedEvent) {
        		removeData(data);
        	}
        	
        	if(treeModel.getChildCount(treeModel.getRoot()) >= 1){
        		cardLayout.last(cardParent);
        		tree.repaint();
        	} else {
        		cardLayout.first(cardParent);
        	}
        }        
	}

	public void valueChanged(TreeSelectionEvent e) {
		if(!disableSelectionReporting){
			
			application.getSelectionManager().clearAll(false, this);

			TreePath[] paths = tree.getSelectionPaths();

			Collection<DataItem> selected = new ArrayList<DataItem>();
			
			for (int i=0; i<paths.length; i++) {			

				DefaultMutableTreeNode node = 
					(DefaultMutableTreeNode) paths[i].getLastPathComponent();
				
				selected.add((DataItem) node.getUserObject());
				
			}
			
			application.getSelectionManager().selectMultiple(selected, this);
		}
	}
}
