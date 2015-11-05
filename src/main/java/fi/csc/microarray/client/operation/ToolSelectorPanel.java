package fi.csc.microarray.client.operation;

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.GridLayout;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
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

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.constants.VisualConstants;

/**
 * The panel for two JLists: on the left side, high-level tool category
 * choice list, and on the right side, a list of lower-level tool choices.
 * Selections in the left list will directly affect the options shown in the
 * right one. Selecting an operation on the right side list, then, will
 * result in corresponding details being shown in the parent ToolPanel.
 * 
 * @author Janne Käki, Aleksi Kallio
 *
 */
public class ToolSelectorPanel extends JPanel
								  implements ListSelectionListener {
	
	private final SwingClientApplication application = (SwingClientApplication)Session.getSession().getApplication();
	
	private ToolPanel toolPanel;
	
	private JList categoryList;
	private JList toolList;

	private ToolModule toolModule;
	
	/**
	 * Creates a new ToolSelectorPanel.
	 * 
	 * @param parent The ToolPanel, for communication purposes.
	 */
	public ToolSelectorPanel(ToolPanel parent, ToolModule toolModule) {
		super(new GridLayout(1, 2));
		this.toolPanel = parent;
		this.toolModule = toolModule;
		
        List<ToolCategory> toolCategories;
        toolCategories = Collections.list(Collections.enumeration(toolModule.getVisibleCategories()));
        
		categoryList = new JList();
		categoryList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		categoryList.addListSelectionListener(this);
		categoryList.setCellRenderer(new CategoryListRenderer());
		categoryList.getInsets().right = 1;
		categoryList.setName("categoryList");
		categoryList.setListData(toolCategories.toArray());
		
		
		
		toolList = new JList();
		toolList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		toolList.addListSelectionListener(this);
		toolList.setCellRenderer(new FontSizeFriendlyListRenderer());
		toolList.addMouseListener(new MouseClickListener());
		toolList.getInsets().right = 1;
		toolList.setName("toolList");
		
		JScrollPane categoryListScroller = new JScrollPane(categoryList);		
		JScrollPane toolListScroller = new JScrollPane(toolList);
		
		//Remove useless borders
		categoryListScroller.setBorder(BorderFactory.createMatteBorder(0, 0, 0, 1,
		        VisualConstants.TOOL_LIST_BORDER_COLOR));
		toolListScroller.setBorder(BorderFactory.createEmptyBorder(0, 0, 0, 0));
		
		this.add(categoryListScroller);
		this.add(toolListScroller);
	}
	
	public Vector<Component> getFocusComponents(){
		Vector<Component> order = new Vector<Component>();
		order.add(categoryList);
		order.add(toolList);		
		return order;
	}
	
	/**
	 * Deselect tool.
	 */
	public void deselectTool() {
	    categoryList.clearSelection();
	    toolList.clearSelection();
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
	

	static class CategoryListRenderer extends FontSizeFriendlyListRenderer {

		public Component getListCellRendererComponent(
				JList list, Object value, int index,
				boolean isSelected, boolean cellHasFocus) {
			
			JLabel comp = (JLabel)super.getListCellRendererComponent(
					list, value, index, isSelected, cellHasFocus);
			
			Color circleColor = ((ToolCategory)(list.getModel().getElementAt(index))).getColor();
			if (circleColor == null) {
				circleColor = comp.getBackground();
			}
			comp.setIcon(new ColoredCircleIcon(circleColor));
			
			return this;
		}
	}
	
	public class ToolPopupMenu extends JPopupMenu {
		JMenuItem helpMenuItem;
		OperationDefinition tool;
		
		public ToolPopupMenu(OperationDefinition tool) {
			this.tool = tool;
			this.helpMenuItem = new JMenuItem("Help...");
			
			helpMenuItem.addActionListener(new ActionListener(){
				public void actionPerformed(ActionEvent e) {
					if (e.getSource() == helpMenuItem){
						application.viewHelpFor(ToolPopupMenu.this.tool);
					}
				}
				
			});
			
			this.add(helpMenuItem);
		}
	}
	
	public class MouseClickListener extends MouseAdapter {
		
		public void mousePressed(MouseEvent e) {
			if (e.getButton() == MouseEvent.BUTTON1) {
				if (e.getClickCount() == 2) {
					toolPanel.executeCurrentOperation(true);
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
					int index = toolList.locationToIndex(e.getPoint());

					// Get component by index
					Object clicked = toolList.getModel().getElementAt(index);

					// Check that cell really is on a clicked point
					Rectangle bounds = toolList.getCellBounds(index, index);

					if(bounds != null && clicked != null) {
						boolean validPoint = bounds.contains(e.getPoint());

						if (clicked instanceof OperationDefinition && validPoint){
							OperationDefinition tool = (OperationDefinition)clicked;

							JPopupMenu popup = new ToolPopupMenu(tool);
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
			if (selected instanceof ToolCategory) {
				ToolCategory selectedCategory = (ToolCategory) selected;
				toolList.clearSelection();
				toolPanel.selectTool(null);
				toolList.setListData(selectedCategory.getToolList());

				// select the category again as selectTool(null) above has cleared the selection
				// disable events during this selection
				categoryList.removeListSelectionListener(this);
        		categoryList.setSelectedValue(selectedCategory, true);
        		categoryList.addListSelectionListener(this);

 			}
		} else if (source == toolList) {
			Object selected = toolList.getSelectedValue();
			if (selected instanceof ExecutionItem) {
				ExecutionItem selectedTool = (ExecutionItem) selected;
				toolPanel.selectTool(selectedTool);
			}
		}
	}

	public String getModuleName() {
		return toolModule.getModuleName();
	}

	/**
	 * 
	 * 
	 * @param category
	 */
	public void selectCategory(ToolCategory category) {
		categoryList.setSelectedValue(category, true);
	}

	public void selectTool(OperationDefinition tool) {
		toolList.setSelectedValue(tool, true);
	}
	
}
