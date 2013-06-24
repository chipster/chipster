package fi.csc.chipster.web.model;

import com.vaadin.data.Container;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.event.DataBoundTransferable;
import com.vaadin.event.dd.DragAndDropEvent;
import com.vaadin.event.dd.DropHandler;
import com.vaadin.event.dd.acceptcriteria.AcceptCriterion;
import com.vaadin.event.dd.acceptcriteria.SourceIs;
import com.vaadin.shared.ui.dd.VerticalDropLocation;
import com.vaadin.ui.AbstractSelect;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.ComboBox;
import com.vaadin.ui.Component;
import com.vaadin.ui.Field;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Table;
import com.vaadin.ui.Table.TableDragMode;
import com.vaadin.ui.TableFieldFactory;
import com.vaadin.ui.TextField;
import com.vaadin.ui.themes.BaseTheme;

import fi.csc.chipster.web.tooledit.Icon;
import fi.csc.chipster.web.tooledit.ToolEditor;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLSyntax.ParameterType;


public class Parameter extends BasicModel{

	private static final long serialVersionUID = 56345103682438254L;
	
	private static final String COLUMN_ID = "ID";
	private static final String COLUMN_NAME = "Name";
	private static final String COLUMN_ACTION = "Action";
	
	private Label lbType;
	private Label lbDefaultValue;
	private Label lbMaxValue;
	private Label lbMinValue;
	
	private ComboBox type;
	private TextField defaultValue;
	private TextField maxValue;
	private TextField minValue;
	private Table typeTable;
	private HorizontalLayout hLayoutTypetable;
	
	public Parameter(ToolEditor root) {
		this.root = root;
		createUI();
	}
	
	@Override
	public Parameter createUI() {
		
//		grid.addComponent(new Label("Parameter"), 0, 0);
		generateHeader();
		initElements();
		addRow(lbType, type);
		addRow(lbMaxValue, maxValue);
		addRow(lbMinValue, minValue);
		addRow(lbDefaultValue, defaultValue);
		addRow(lbOptional, optional);
		addRow(lbDescription, description);
		generateFooter();
		return this;
	}
	
	private void initElements() {
		
		lbType = new Label("Type:");
		lbOptional = new Label("Optional:");
		lbDefaultValue = new Label("Default:");
		lbMaxValue = new Label("Maximum value:");
		lbMinValue = new Label("Minimum value:");
		lbId.setValue("User defined parameter");
		
		type = new ComboBox();
		type.setImmediate(true);
		type.setWidth(WIDTH);
		for(ParameterType parameterType : ParameterType.values()) {
			type.addItem(parameterType);
		}
		type.setNullSelectionAllowed(false);
		type.select(type.getItemIds().iterator().next());
		type.addValueChangeListener(new ValueChangeListener() {
			
			/**
			 * 
			 */
			private static final long serialVersionUID = 3675044736996182900L;

			@Override
			public void valueChange(ValueChangeEvent event) {
				if(type.getValue().equals(ParameterType.ENUM)) {
					addTypeTable();
				} else {
					removeTypeTable();
				}
				String text = getValue(getNameValue()) + " " + getValue(getTypeValue());
				setTitleDescriptionValue(text.trim());
				root.getToolEditorUI().getTreeToolEditor().setItemCaption(Parameter.this, Parameter.this.toString());
			}
		});
		defaultValue = new TextField();
		defaultValue.setWidth(WIDTH);
		maxValue = new TextField();
		maxValue.setWidth(WIDTH);
		minValue = new TextField();
		minValue.setWidth(WIDTH);
		initTypeTable();
	}
	
	public Parameter createUIWithData(SADLDescription.Parameter parameter) {
//		createUI();
		fillWithData(parameter);
		return this;
	}
	
	public void fillWithData(SADLDescription.Parameter parameter) {
		parameter.getType();
		type.select(parameter.getType());
		
		if(parameter.getSelectionOptions() != null) {
//			initTypeTable();
			fillTypeTableWithData(parameter.getSelectionOptions());
			addTypeTable();
		}
		id.setValue(getValue(parameter.getName().getID()));
		name.setValue(getValue(parameter.getName().getDisplayName()));
		description.setValue(getValue(parameter.getComment()));
		optional.setValue(parameter.isOptional());
		if(parameter.getDefaultValues().length >= 1)
			defaultValue.setValue(parameter.getDefaultValues()[0]);
		maxValue.setValue(getValue(parameter.getTo()));
		minValue.setValue(getValue(parameter.getFrom()));
	}
	
	private void initTypeTable() {
		hLayoutTypetable = new HorizontalLayout();
		hLayoutTypetable.setSpacing(true);
		typeTable = new Table();
		typeTable.addContainerProperty(COLUMN_ID, String.class, "");
		typeTable.addContainerProperty(COLUMN_NAME, String.class, "");
		typeTable.addContainerProperty(COLUMN_ACTION, Button.class, null);
		typeTable.setPageLength(1);
		typeTable.setEditable(true);
		typeTable.setImmediate(true);
		typeTable.setTableFieldFactory(new TableFieldFactory() {
			private static final long serialVersionUID = -5275144703109311990L;

			@SuppressWarnings("unchecked")
			@Override
			public Field<?> createField(Container container, final Object itemId,
					final Object propertyId, Component uiContext) {
					
				final TextField txtField = new TextField();
				if(propertyId.equals(COLUMN_ID)) {
					txtField.setValue(((Name) itemId).getID());
					typeTable.getContainerProperty(itemId, COLUMN_ID).setValue(getValue(((Name) itemId).getID()));
				} else if(propertyId.equals(COLUMN_NAME)) {
					txtField.setValue(((Name) itemId).getDisplayName());
					typeTable.getContainerProperty(itemId, COLUMN_NAME).setValue(getValue(((Name) itemId).getDisplayName()));
				}
				
				txtField.addValueChangeListener(new ValueChangeListener() {
					
					/**
					 * 
					 */
					private static final long serialVersionUID = 7785711696111601557L;

					@Override
					public void valueChange(ValueChangeEvent event) {
						if(propertyId.equals(COLUMN_ID)) {
							((Name) itemId).setID(txtField.getValue());
						} else if(propertyId.equals(COLUMN_NAME)) {
							((Name) itemId).setDisplayName(txtField.getValue());
						}
					}
				});
				
				if(typeTable.getContainerDataSource().size() > 5) {
					typeTable.setPageLength(5);
				} else {
					typeTable.setPageLength(typeTable.getContainerDataSource().size());
				}
				
				return txtField;
			}
			
		});
		typeTable.setDragMode(TableDragMode.ROW);
		typeTable.setDropHandler(new DropHandler() {
			private static final long serialVersionUID = 1L;

			@Override
			public AcceptCriterion getAcceptCriterion() {
				return new SourceIs(typeTable);
			}
			
			@SuppressWarnings("unchecked")
			@Override
			public void drop(DragAndDropEvent event) {
				DataBoundTransferable transferable = (DataBoundTransferable) event.getTransferable();
		        Object sourceRow = transferable.getItemId();

		        AbstractSelect.AbstractSelectTargetDetails dropData =
		                ((AbstractSelect.AbstractSelectTargetDetails) event.getTargetDetails());
		        Object targetRow = dropData.getItemIdOver();

		        if ((sourceRow == targetRow) || (targetRow == null)) {
		            return;
		        }
		        typeTable.removeItem(sourceRow);
		        if (dropData.getDropLocation() == VerticalDropLocation.BOTTOM) {
		            typeTable.addItemAfter(targetRow, sourceRow);
		        }
		        else {
		            Object rowAbove = typeTable.prevItemId(targetRow);
		            typeTable.addItemAfter(rowAbove, sourceRow);
		        }
		        typeTable.getContainerProperty(sourceRow, COLUMN_ACTION).setValue(getDeleteButton(sourceRow));
			}
		});
		
		hLayoutTypetable.addComponent(typeTable);
		Button btAddNewRow = new Button("New Row");
		btAddNewRow.setIcon(Icon.getResource(Icon.getAddButtonIconPath()));
		btAddNewRow.addClickListener(new ClickListener() {
			private static final long serialVersionUID = 2039910564036291751L;

			@SuppressWarnings("unchecked")
			@Override
			public void buttonClick(ClickEvent event) {
				Name name = Name.createEmptyName();
				typeTable.addItem(name);
				typeTable.getContainerProperty(name, COLUMN_ACTION).setValue(getDeleteButton(name));
				typeTable.setValue(name);
				typeTable.setCurrentPageFirstItemId(name);
			}
		});
		hLayoutTypetable.addComponent(btAddNewRow);
	}
	
	private void fillTypeTableWithData(final Name[] types) {
		for(int i = 0; i < types.length; i++) {
			typeTable.addItem(new Object[] {types[i].getID(), types[i].getDisplayName(), getDeleteButton(types[i])}, types[i]);
		}
//		if(typeTable.getContainerDataSource().size() > 5) {
//			typeTable.setPageLength(5);
//		} else {
//			typeTable.setPageLength(typeTable.getContainerDataSource().size());
//		}
	}
	
	private void addTypeTable() {
		if(this.getComponentArea(hLayoutTypetable) == null) {
			int i = 0;
			while(!this.getComponent(1, i).equals(type)) {
				i++;
			}
			this.insertRow(i+1);
			this.addComponent(hLayoutTypetable, 1, i+1);
		}
	}
	
	private void removeTypeTable() {
		this.removeComponent(hLayoutTypetable);
	}
	
	public SADLDescription.Parameter getSadlParameter() {
		SADLDescription.Parameter parameter = new SADLDescription.Parameter(Name.createName(id.getValue(), name.getValue()), (ParameterType) type.getValue(), 
				getEnumList(), getValueOrNull(minValue.getValue()), 
				getValueOrNull(maxValue.getValue()), getValueOrNull(defaultValue.getValue()), getValueOrNull(description.getValue()));
		parameter.setOptional(optional.getValue());
		return parameter;
	}
	
	private Name[] getEnumList() {
		if(typeTable == null || typeTable.getItemIds().isEmpty())
			return null;
		Name[] names = typeTable.getItemIds().toArray(new Name[0]);
		return names;
	}
	
	private Button getDeleteButton(final Object itemId) {
		Button btDelete = new Button();
		btDelete.setIcon(Icon.getResource(Icon.getDeleteButtonIconPath()));
		btDelete.setStyleName(BaseTheme.BUTTON_LINK);
		btDelete.addClickListener(new ClickListener() {
			private static final long serialVersionUID = -3695725710938486562L;

			@Override
			public void buttonClick(ClickEvent event) {
				typeTable.removeItem(itemId);
				
			}
		});
		return btDelete;
	}

	@Override
	protected void generateHeader() {
		lbTitle.setValue(getBoldText("Parameter"));
	}

	@Override
	protected String getType() {
		// TODO Auto-generated method stub
		return type.getValue().toString();
	}
	
	@Override
	public String toString() {
		return getValue(name.getValue()) + " " + getValue(type.getValue().toString()) + " " + (optional.getValue() ? " OPTIONAL" : "");
	}
}
