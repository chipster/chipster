package fi.csc.chipster.web.tooledit;

import java.util.Collection;

import com.vaadin.data.util.HierarchicalContainer;
import com.vaadin.event.DataBoundTransferable;
import com.vaadin.event.ItemClickEvent;
import com.vaadin.event.ItemClickEvent.ItemClickListener;
import com.vaadin.event.Transferable;
import com.vaadin.event.dd.DragAndDropEvent;
import com.vaadin.event.dd.DropHandler;
import com.vaadin.event.dd.acceptcriteria.AcceptCriterion;
import com.vaadin.event.dd.acceptcriteria.Not;
import com.vaadin.event.dd.acceptcriteria.Or;
import com.vaadin.server.ThemeResource;
import com.vaadin.shared.ui.dd.VerticalDropLocation;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.Component;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Tree;
import com.vaadin.ui.TreeTable;
import com.vaadin.ui.themes.BaseTheme;

import fi.csc.chipster.web.model.BasicModel;
import fi.csc.chipster.web.model.Input;
import fi.csc.chipster.web.model.Output;
import fi.csc.chipster.web.model.Parameter;
import fi.csc.chipster.web.model.Tool;

/**
 * Summary tree for adding/deleting/changing positions(drag & drop) inputs, outputs, parameters
 * @author Gintare Pacauskaite
 *
 */
public class TreeToolEditor extends TreeTable implements ItemClickListener{
	
	private static final long serialVersionUID = -8713542953586370559L;
	
	private static final String ACTIONS = "Actions";
	
	private ToolEditorUI root;
	
	private static final String TOOL = "Tool";
	private static final String INPUTS = "Inputs";
	private static final String OUTPUTS = "Outputs";
	private static final String PARAMETERS = "Parameters";
	
	private Tool tool;
	
	
	public TreeToolEditor(ToolEditorUI root) {
		this.root = root;
		tool = new Tool(root.getToolEditor());
		initTree();
	}
	
	private void initTree() {
		this.setSizeFull();
		this.setImmediate(true);
		this.setSelectable(true);
		this.addItemClickListener(this);
		this.addContainerProperty("", String.class, null);
		this.addContainerProperty(ACTIONS, HorizontalLayout.class, null	);
		this.addItem(new Object[] {TOOL, getActionLayout(false, false, TOOL)}, TOOL);
		this.setItemDescriptionGenerator(new ItemDescriptionGenerator() {
			private static final long serialVersionUID = -1913286695570843896L;

			@Override
			public String generateDescription(Component source, Object itemId,
					Object propertyId) {
				String description = "Show all ";
				if(itemId.equals(TOOL))
					description += "elements";
				else if (itemId.equals(INPUTS))
					description += "inputs";
				else if (itemId.equals(OUTPUTS))
					description += "outputs";
				else if (itemId.equals(PARAMETERS))
					description += "parameters";
				else
					return null;
				return description;
			}
		});
		this.setCollapsed(TOOL, false);
		String[] rootToolElements = new String[] {INPUTS, OUTPUTS, PARAMETERS};
		for(String element : rootToolElements) {
			this.addItem(new Object[] {element, getActionLayout(true, false, element)}, element);
			this.setParent(element, TOOL);
			this.setCollapsed(element, false);
		}
		this.setDragMode(TableDragMode.ROW);
		this.setDropHandler(new DropHandler() {
			private static final long serialVersionUID = -4415321436294383112L;

			@Override
			public AcceptCriterion getAcceptCriterion() {
				return new Or(Tree.TargetItemAllowsChildren.get(), new Not(VerticalLocationIs.MIDDLE));
			}
			
			@Override
			public void drop(DragAndDropEvent event) {
	            final Transferable t = event.getTransferable();
	            if (t.getSourceComponent() != TreeToolEditor.this
	                    || !(t instanceof DataBoundTransferable)) {
	                return;
	            }

	            final AbstractSelectTargetDetails dropData = ((AbstractSelectTargetDetails) event.getTargetDetails());

	            final Object sourceItemId = ((DataBoundTransferable) t).getItemId();
	            final Object targetItemId = dropData.getItemIdOver();
	            final VerticalDropLocation location = dropData.getDropLocation();

	            moveNode(sourceItemId, targetItemId, location);
			}
		});

	}
	
	/**
	 * drag and drop 
	 * @param sourceItemId what was dragged
	 * @param targetItemId on what was dropped
	 * @param location bottom, top, middle
	 */
	private void moveNode(final Object sourceItemId,
            final Object targetItemId, final VerticalDropLocation location) {
        final HierarchicalContainer container = (HierarchicalContainer) this.getContainerDataSource();
        if(sourceItemId instanceof String)
        	return;
        if((sourceItemId instanceof Input) && (!(targetItemId instanceof Input) && !targetItemId.equals(INPUTS)))
        	return;
        if((sourceItemId instanceof Output) && (!(targetItemId instanceof Output) && !targetItemId.equals(OUTPUTS)))
        	return;
        if((sourceItemId instanceof Parameter) && (!(targetItemId instanceof Parameter) && !targetItemId.equals(PARAMETERS)))
        	return;
        
        if (location == VerticalDropLocation.MIDDLE) {
            if (container.setParent(sourceItemId, targetItemId)
                    && container.hasChildren(targetItemId)) {
                container.moveAfterSibling(sourceItemId, null);
            }
        } else if (location == VerticalDropLocation.TOP) {
            final Object parentId = container.getParent(targetItemId);
            if(parentId.equals(TOOL)) {
            	Object lastChild = container.getChildren(targetItemId).toArray()[container.getChildren(targetItemId).size()-1];
            	if(!lastChild.equals(sourceItemId)) {
            		container.moveAfterSibling(sourceItemId, lastChild);
            	}
            } else if (container.setParent(sourceItemId, parentId)) {
                // reorder only the two items, moving source above target
                container.moveAfterSibling(sourceItemId, targetItemId);
                container.moveAfterSibling(targetItemId, sourceItemId);
            } 
        } else if (location == VerticalDropLocation.BOTTOM) {
            final Object parentId = container.getParent(targetItemId);
            if (parentId.equals(TOOL)) {
            	Object firstChild = container.getChildren(targetItemId).iterator().next();
            	if(!firstChild.equals(sourceItemId)) {
            		container.moveAfterSibling(sourceItemId, firstChild);
            		container.moveAfterSibling(firstChild, sourceItemId);
            	}
            } else if (container.setParent(sourceItemId, parentId)) {
                container.moveAfterSibling(sourceItemId, targetItemId);
            }
        }
    }
	
	private HorizontalLayout getActionLayout(boolean needAddButton, final Object itemId) {
		return getActionLayout(needAddButton, true, itemId);
	}
	
	private HorizontalLayout getActionLayout(boolean needAddButton, boolean needDeleteButton, final Object itemId) {
		HorizontalLayout hLayout = new HorizontalLayout();
		hLayout.setSpacing(true);
		if(needAddButton) {
			Button btAdd = new Button();
			hLayout.addComponent(btAdd);
			btAdd.setIcon(new ThemeResource("images/edit-add-2.png"));
			btAdd.setStyleName(BaseTheme.BUTTON_LINK);
			btAdd.setDescription("Add " + (itemId.equals(INPUTS) ? "Input" : (itemId.equals(OUTPUTS) ? "Output" : "Parameter")));
			
			btAdd.addClickListener(new ClickListener() {
				private static final long serialVersionUID = 1L;

				@Override
				public void buttonClick(ClickEvent event) {
					BasicModel model = getModelByParent(itemId);
					addElement(model, itemId);
					TreeToolEditor.this.select(model);
					root.getToolEditor().removeAllComponents();
					root.getToolEditor().addComponent(model);
				}
			});
		}
		if(needDeleteButton) {
			Button btDelete = new Button();
			btDelete.setIcon(new ThemeResource("images/close.png"));
			btDelete.setStyleName(BaseTheme.BUTTON_LINK);
			String description = "Delete ";
			description += (itemId instanceof Input ? "Input" : (itemId instanceof Output ? "Output" : "Parameter"));
			btDelete.setDescription(description);
			hLayout.addComponent(btDelete, 0);
			btDelete.addClickListener(new ClickListener() {
				private static final long serialVersionUID = 1L;
	
				@Override
				public void buttonClick(ClickEvent event) {
					if(itemId.equals(TOOL)) {
						removeAllChildren();
					} else if(itemId.equals(INPUTS)) {
						removeChildren(INPUTS);
					} else if(itemId.equals(OUTPUTS)) {
						removeChildren(OUTPUTS);
					} else if(itemId.equals(PARAMETERS)) {
						removeChildren(PARAMETERS);
					} else {
						root.getToolEditor().removeComponent((BasicModel) itemId);
						TreeToolEditor.this.removeItem(itemId);
					}
				}
			});
		}

		return hLayout;
	}
	
	private BasicModel getModelByParent(Object itemId) {
		return (itemId.equals(INPUTS)) ? new Input(root.getToolEditor()) : ((itemId.equals(OUTPUTS)) ? new Output(root.getToolEditor()) : new Parameter(root.getToolEditor()));
	}

	@Override
	public void itemClick(ItemClickEvent event) {
		Object itemId = event.getItemId();
		if(itemId instanceof Input || itemId instanceof Output || itemId instanceof Parameter) {
			root.getToolEditor().removeAllComponents();
			root.getToolEditor().addComponent((BasicModel)itemId);
		} else {
			root.getToolEditor().removeAllComponents();
			if(itemId.equals(TOOL)) {
				root.getToolEditor().addComponent(tool);
				addAllChildrenToToolEditor(INPUTS);
				addAllChildrenToToolEditor(OUTPUTS);
				addAllChildrenToToolEditor(PARAMETERS);
			} else if(itemId.equals(INPUTS)) {
				addAllChildrenToToolEditor(INPUTS);
			} else if(itemId.equals(OUTPUTS)) {
				addAllChildrenToToolEditor(OUTPUTS);
			} else if(itemId.equals(PARAMETERS)) {
				addAllChildrenToToolEditor(PARAMETERS);
			}
		}
	}
	
	@SuppressWarnings("unchecked")
	private void addAllChildrenToToolEditor(Object parent) {
		Collection<BasicModel> elements = (Collection<BasicModel>) this.getChildren(parent);
		if(elements != null)
			for(BasicModel element : elements) {
				root.getToolEditor().addComponent(element);
			}
	}
	
	@SuppressWarnings("unchecked")
	public void setItemCaption(Object itemId, String caption) {
		try {
			this.getContainerProperty(itemId, "").setValue(caption);
	
		} catch(NullPointerException e) {
			
		}
	}
	
	public String getItemCaption(Object itemId, String caption) {
		return caption;
	}
	
	private void removeChildren(Object parent) {
		@SuppressWarnings("unchecked")
		Collection<BasicModel> elements = (Collection<BasicModel>) this.getChildren(parent);
		if(elements == null || elements.isEmpty())
			return;
		for(BasicModel element : elements.toArray(new BasicModel[0])) {
			root.getToolEditor().removeComponent(element);
			this.removeItem(element);
		}
	}
	
	@SuppressWarnings("unchecked")
	public Collection<Input> getInputs() {
		return (Collection<Input>) this.getChildren(INPUTS);
	}
	
	@SuppressWarnings("unchecked")
	public Collection<Output> getOutputs() {
		return (Collection<Output>) this.getChildren(OUTPUTS);
	}
	
	@SuppressWarnings("unchecked")
	public Collection<Parameter> getParameters() {
		return (Collection<Parameter>) this.getChildren(PARAMETERS);
	}
	
	public void removeAllChildren() {
		tool = new Tool(root.getToolEditor());
		removeChildren(INPUTS);
		removeChildren(OUTPUTS);
		removeChildren(PARAMETERS);
	}
	
	public void addInput(Input input) {
		addElement(input, INPUTS);
	}
	
	public void addOutput(Output output) {
		addElement(output, OUTPUTS);
	}
	
	public void addParameter(Parameter parameter) {
		addElement(parameter, PARAMETERS);
	}
	
	public void addElement(BasicModel model, Object itemId) {
		this.addItem(new Object[] {(Object) getItemCaption(model, model.toString()), getActionLayout(false, model)},  model);
		model.setTitleDescriptionValue(model.toString());
		this.setChildrenAllowed(model, false);
		this.setParent(model, itemId);
	}
	
	public void addTool(Tool tool) {
		this.tool = tool;
	}
	
	public Tool getTool() {
		return tool;
	}
	
	@SuppressWarnings("unchecked")
	public void updateToolTitle(String text) {
		getContainerProperty(TOOL, "").setValue((text.isEmpty() ? "Tool" : text));
	}
}
