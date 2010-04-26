package fi.csc.microarray.client.operation;

import java.awt.Component;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.DefaultListCellRenderer;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;
import javax.swing.UIManager;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;

/**
 * The panel for two JLists: on the left side, high-level operation category
 * choice list, and on the right side, a list of lower-level operation choices.
 * Selections in the left list will directly affect the options shown in the
 * right one. Selecting an operation on the right side list, then, will
 * result in corresponding details being shown in the parent OperationPanel.
 * 
 * @author Janne KÃ¤ki
 *
 */
@SuppressWarnings("serial")
public class OperationChoicePanel extends JPanel
								  implements ListSelectionListener {
	// Logger for this class
	private static final Logger logger = Logger
			.getLogger(OperationChoicePanel.class);
	
	private final ClientApplication application = Session.getSession().getApplication();
	
	private OperationPanel parent;
	
	private JList categoryList;
	private JList operationList;
	
	private OperationCategory selectedCategory;
	private ExecutionItem selectedOperation;
	
	/**
	 * Creates a new OperationChoicePanel.
	 * 
	 * @param parent The OperationPanel, for communication purposes.
	 */
	public OperationChoicePanel(OperationPanel parent,
	       Collection<OperationCategory> operationCategoryCollection) {
		super(new GridLayout(1, 2));
		this.parent = parent;

        List<OperationCategory> operationCategories;
        operationCategories = Collections.list(Collections.enumeration(operationCategoryCollection));

		Object[] categories;
		if (operationCategories != null) {
			categories = new Object[operationCategories.size()];
			for (int i = 0; i < operationCategories.size(); i++) {
				categories[i] = operationCategories.get(i);
			}
		} else {
			categories = new Object[0];
		}
		
		categoryList = new JList(categories);
		categoryList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		categoryList.addListSelectionListener(this);
		categoryList.setCellRenderer(new CategoryListRenderer());
		//categoryList.setPreferredSize(new Dimension(130, 0));
		categoryList.getInsets().right = 1;
		categoryList.setName("categoryList");
		
		operationList = new JList();
		operationList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		operationList.addListSelectionListener(this);
		operationList.setCellRenderer(new FontSizeFriendlyListRenderer());
		operationList.addMouseListener(new MouseClickListener());
		//operationList.setPreferredSize(new Dimension(200, 0));
		operationList.getInsets().right = 1;
		operationList.setName("operationList");
		
		JScrollPane categoryListScroller = new JScrollPane(categoryList);		
		JScrollPane operationListScroller = new JScrollPane(operationList);
		
		//Remove useless borders
		categoryListScroller.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
		operationListScroller.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
		
		this.add(categoryListScroller);
		this.add(operationListScroller);
	}
	
	public Vector<Component> getFocusComponents(){
		Vector<Component> order = new Vector<Component>();
		order.add(categoryList);
		order.add(operationList);		
		return order;
	}
	
	/**
	 * Deselect operation.
	 */
	public void deselectOperation() {
	    categoryList.clearSelection();
	    operationList.clearSelection();
	    parent.selectOperation(null);
	}
	
	static class FontSizeFriendlyListRenderer extends DefaultListCellRenderer {
		public Component getListCellRendererComponent(
				JList list, Object value, int index,
				boolean isSelected, boolean cellHasFocus) {
			
			JLabel comp = (JLabel)super.getListCellRendererComponent(
					list, value, index, isSelected, cellHasFocus);

			Font font = UIManager.getFont("Label.font");
			comp.setFont(font);			
			
			return this;
		}
	}
	

	class CategoryListRenderer extends FontSizeFriendlyListRenderer {

		public Component getListCellRendererComponent(
				JList list, Object value, int index,
				boolean isSelected, boolean cellHasFocus) {
			
			JLabel comp = (JLabel)super.getListCellRendererComponent(
					list, value, index, isSelected, cellHasFocus);
			
			comp.setIcon(new ColoredCircleIcon(((OperationCategory)(list.getModel().getElementAt(index))).getColor()));
			
			return this;
		}
	}
	
	public class OperationPopupMenu extends JPopupMenu {
		JMenuItem helpMenuItem;
		OperationDefinition operationDefinition;
		
		public OperationPopupMenu(OperationDefinition operationDefinition) {
			this.operationDefinition = operationDefinition;
			this.helpMenuItem = new JMenuItem("Help...");
			
			helpMenuItem.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (e.getSource() == helpMenuItem){
						application.viewHelpFor(OperationPopupMenu.this.operationDefinition);
					}
				}
				
			});
			
			this.add(helpMenuItem);
		}
	}
	
	public class MouseClickListener extends MouseAdapter {
		
		public void mousePressed(MouseEvent e) {
			if (e.getButton() == MouseEvent.BUTTON1) {
				if (e.getClickCount() > 1) {
					Object selected = operationList.getSelectedValue();
					if (selected instanceof OperationDefinition &&
							application.getSelectionManager().getSelectedDataBeans().size() > 0) {
						
						OperationDefinition operationDefinition = (OperationDefinition)selected;
						if (!operationDefinition.evaluateSuitabilityFor(application.getSelectionManager().getSelectedDataBeans()).isImpossible()) {
							application.executeOperation(operationDefinition, null);
						}
					}
				}
			}
			maybeShowPopup(e);        
		}	
		@Override
		public void mouseReleased(MouseEvent e) {
			maybeShowPopup(e);        
		}

		private void maybeShowPopup(MouseEvent e) {
			if (e.isPopupTrigger()) {
//				Help popup menu
				if (e.getButton() == MouseEvent.BUTTON3) {
					// Get index of the clicked component
					int index = operationList.locationToIndex(e.getPoint());

					// Get component by index
					Object clicked = operationList.getModel().getElementAt(index);

					// Check that cell really is on a clicked point
					Rectangle bounds = operationList.getCellBounds(index, index);

					if(bounds != null && clicked != null) {
						boolean validPoint = bounds.contains(e.getPoint());

						if (clicked instanceof OperationDefinition && validPoint){
							OperationDefinition operationDefinition = (OperationDefinition)clicked;

							JPopupMenu popup = new OperationPopupMenu(operationDefinition);
							popup.show(e.getComponent(), e.getX(), e.getY());
						}
					}
				}
			}
		}
	}

	
	
	/**
	 * A method which allows the panel to interact with the user's selections.
	 * Defined by the ListSelectionlistener interface.
	 */
	public void valueChanged(ListSelectionEvent e) {
		Object source = e.getSource();
		if (source == categoryList) {
			Object selected = categoryList.getSelectedValue();
			if (selected instanceof OperationCategory) {
				selectedCategory = (OperationCategory) selected;
				selectedOperation = null;
				logger.debug("selected category has " + selectedCategory.getOperationList().size() + " operations");
				operationList.setListData(selectedCategory.getOperationList());
				parent.enableAction(false);
			}
			parent.selectOperation(null);
		} else if (source == operationList) {
			Object selected = operationList.getSelectedValue();
			if (selected instanceof ExecutionItem) {
				selectedOperation = (ExecutionItem) selected;
				parent.selectOperation(selectedOperation);
			}
		}
	}
}
