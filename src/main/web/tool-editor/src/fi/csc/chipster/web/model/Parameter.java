package fi.csc.chipster.web.model;

import java.util.ArrayList;
import java.util.Collection;

import com.vaadin.data.Container;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.data.Validator.InvalidValueException;
import com.vaadin.data.Validator;
import com.vaadin.event.DataBoundTransferable;
import com.vaadin.event.dd.DragAndDropEvent;
import com.vaadin.event.dd.DropHandler;
import com.vaadin.event.dd.acceptcriteria.AcceptCriterion;
import com.vaadin.event.dd.acceptcriteria.SourceIs;
import com.vaadin.server.Page;
import com.vaadin.server.ThemeResource;
import com.vaadin.shared.ui.dd.VerticalDropLocation;
import com.vaadin.shared.ui.slider.SliderOrientation;
import com.vaadin.ui.AbstractSelect;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.ComboBox;
import com.vaadin.ui.Component;
import com.vaadin.ui.Field;
import com.vaadin.ui.Label;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;
import com.vaadin.ui.Slider;
import com.vaadin.ui.Table;
import com.vaadin.ui.Table.TableDragMode;
import com.vaadin.ui.TableFieldFactory;
import com.vaadin.ui.TextField;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.themes.BaseTheme;

import fi.csc.chipster.web.tooledit.Icon;
import fi.csc.chipster.web.tooledit.ToolEditor;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

/**
 * Tool editor parameter model
 * @author Gintare Pacauskaite
 *
 */
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
	//for enum type
	private ComboBox defaultValue2;
	// for percent type
	private Slider defaultSlider;
	private TextField maxValue;
	private TextField minValue;
	private Table typeTable;
	private VerticalLayout vLayoutTypetable;
	
	public Parameter(ToolEditor root) {
		this.root = root;
		createUI();
	}
	
	@Override
	public Parameter createUI() {
		
		generateHeader();
		initElements();
		addRow(lbType, type);
		addRow(lbMaxValue, maxValue);
		addRow(lbMinValue, minValue);
		addRow(lbDefaultValue, defaultValue);
		addRow(lbOptional, optional);
		addRow(lbDescription, description);
		return this;
	}
	
	private void initElements() {
		
		lbType = new Label("Type:");
		lbOptional = new Label("Optional:");
		lbDefaultValue = new Label("Default:");
		lbMaxValue = new Label("Maximum value:");
		lbMinValue = new Label("Minimum value:");
		lbId.setValue("Id:");
		
		type = new ComboBox();
		type.setImmediate(true);
		type.setWidth(WIDTH);
		for(ParameterType parameterType : ParameterType.values()) {
			type.addItem(parameterType);
		}
		type.setNullSelectionAllowed(false);
		type.select(type.getItemIds().iterator().next());
		type.addValueChangeListener(new ValueChangeListener() {
			private static final long serialVersionUID = 3675044736996182900L;

			@Override
			public void valueChange(ValueChangeEvent event) {
				defaultValue.setValue("");
				minValue.setValue("");
				maxValue.setValue("");
				if(type.getValue().equals(ParameterType.ENUM)) {
					addTypeTable();
					changeDefautlToDropDown();
					fillDefaultValueDropdownWithTypeData();
					removeMinMaxFields();
				} else {
					removeTypeTable();
					changeDefaultToTextField();
					AddMinMaxFields();
					if(type.getValue().equals(ParameterType.STRING) || type.getValue().equals(ParameterType.COLUMN_SEL)
							|| type.getValue().equals(ParameterType.INPUT_SEL) || type.getValue().equals(ParameterType.METACOLUMN_SEL)
							|| type.getValue().equals(ParameterType.PERCENT)) {
						removeMinMaxFields();
					}
					if(type.getValue().equals(ParameterType.PERCENT)) {
						ChangeDefaultToSlider();
					}
				}
				String text = getValue(getNameValue()) + " " + getValue(getTypeValue());
				setTitleDescriptionValue(text.trim());
				root.getToolEditorUI().getTreeToolEditor().setItemCaption(Parameter.this, Parameter.this.toString());
			}
		});
		
		defaultValue = new TextField();
		defaultValue.setDescription("Default value of parameter. Can be empty");
		defaultValue.setImmediate(true);
		defaultValue.addValidator(new Validator() {
			private static final long serialVersionUID = 1739911601024686605L;

			@Override
			public void validate(Object value) throws InvalidValueException {
				if(value.toString().isEmpty())
					return;
				validateNumberForNumberType(value);
				double dvalue = Double.parseDouble(value.toString());
				if(!maxValue.getValue().isEmpty()) {
					try {
						double max = Double.parseDouble(maxValue.getValue());
						if(dvalue > max) {
							String text = "Value must be smaller than " + max;
							new Notification(text, Type.WARNING_MESSAGE).show(Page.getCurrent());
							throw new InvalidValueException(text);
						}
					} catch(NumberFormatException e) {}
				}
				if(!minValue.getValue().isEmpty()) {
					try {
						double min = Double.parseDouble(minValue.getValue());
						if(dvalue < min) {
							String text = "Value must be greater than " + min;
							new Notification(text, Type.WARNING_MESSAGE).show(Page.getCurrent());
							throw new InvalidValueException(text);
						}
					} catch(NumberFormatException e) {}
				}
			}
		});
		defaultValue.setWidth(WIDTH);
		defaultValue2 = new ComboBox();
		defaultValue2.setImmediate(true);
		defaultValue2.setWidth(WIDTH);
		defaultSlider = new Slider();
		defaultSlider.setWidth(WIDTH);
		defaultSlider.setImmediate(true);
		defaultSlider.setOrientation(SliderOrientation.HORIZONTAL);
		
		maxValue = new TextField();
		maxValue.setDescription("Maximum value of parameter. Can be empty");
		maxValue.setImmediate(true);
		maxValue.addValidator(new Validator() {
			private static final long serialVersionUID = 1739911601024686605L;

			@Override
			public void validate(Object value) throws InvalidValueException {
				validateNumberForNumberType(value);
				if(!minValue.getValue().isEmpty()) {
					try {
						double min = Double.parseDouble(minValue.getValue());
						double max = Double.parseDouble(value.toString());
						if(max < min) {
							throw new InvalidValueException("Maximum value must be greater than minimum value");
						}
					} catch (NumberFormatException e) {}
				}
			}
		});
		maxValue.setWidth(WIDTH);
		minValue = new TextField();
		minValue.setWidth(WIDTH);
		minValue.setImmediate(true);
		minValue.setDescription("Minimum value of parameter. Can be empty");
		minValue.addValidator(new Validator() {
			private static final long serialVersionUID = -448212312378806070L;

			@Override
			public void validate(Object value) throws InvalidValueException {
				validateNumberForNumberType(value);
				if(!maxValue.getValue().isEmpty()) {
					try {
						double max = Double.parseDouble(maxValue.getValue());
						double min = Double.parseDouble(value.toString());
						if(min > max) {
							throw new InvalidValueException("Minimum value must be smaller than maximum value");
						}
					} catch (NumberFormatException e) {}
				}
			}
		});
		initTypeTable();
	}
	
	private void validateNumberForNumberType(Object value) throws InvalidValueException {
		if(value.toString().isEmpty())
			return;
		if(type.getValue().equals(ParameterType.DECIMAL) || type.getValue().equals(ParameterType.INTEGER)) {
			if (!isNumber((String) value)) {
				Notification.show("Value must be a number", Type.ERROR_MESSAGE);
				throw new InvalidValueException("Value must be a number");
			}
		}
	}
	
	private boolean isNumber(String value) {
			try {
				Double.parseDouble(value);
			} catch (NumberFormatException e) {
				return false;
			}
			return true;
	}
	
	/**
	 * Creats parameter from SADL parameter
	 * @param parameter SADL parameter
	 * @return tool editor parameter
	 */
	public Parameter createUIWithData(SADLDescription.Parameter parameter) {
		fillWithData(parameter);
		return this;
	}
	
	/**
	 * Fills up parameter fields with SADL parameter data
	 * @param parameter SADL parameter
	 */
	public void fillWithData(SADLDescription.Parameter parameter) {
		type.select(parameter.getType());
		
		if(parameter.getSelectionOptions() != null) {
			fillTypeTableWithData(parameter.getSelectionOptions());
			addTypeTable();
			fillDefaultValueDropdownWithTypeData();
			
		}
		id.setValue(getValue(parameter.getName().getID()));
		name.setValue(getValue(parameter.getName().getDisplayName()));
		description.setValue(getValue(parameter.getComment()));
		optional.setValue(parameter.isOptional());
		if(parameter.getDefaultValues().length >= 1) {
			String defaultVal = parameter.getDefaultValues()[0];
			if(parameter.getType().equals(ParameterType.ENUM)) {
				Name itemId = getdefaultObject(defaultVal);
				if(name != null) {
					defaultValue2.select(itemId);
				}
			} else if(parameter.getType().equals(ParameterType.PERCENT)) {
				try { 
					defaultSlider.setValue(Double.valueOf(defaultVal));
				} catch(NumberFormatException e) {}
			} else {
				defaultValue.setValue(defaultVal);
			}
		}
		maxValue.setValue(getValue(parameter.getTo()));
		minValue.setValue(getValue(parameter.getFrom()));
	}
	
	@SuppressWarnings("unchecked")
	private Name getdefaultObject(String id) {
		for(Name name : (Collection<Name>) defaultValue2.getItemIds()) {
			if(name.getID().equals(id)) {
				return name;
			}
		}
		return null;
	}
	
	private void initTypeTable() {
		vLayoutTypetable = new VerticalLayout();
		vLayoutTypetable.setSpacing(true);
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
					private static final long serialVersionUID = 7785711696111601557L;

					@Override
					public void valueChange(ValueChangeEvent event) {
						if(propertyId.equals(COLUMN_ID)) {
							((Name) itemId).setID(txtField.getValue());
						} else if(propertyId.equals(COLUMN_NAME)) {
							((Name) itemId).setDisplayName(txtField.getValue());
							defaultValue2.setItemCaption(itemId, txtField.getValue());
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
		
		vLayoutTypetable.addComponent(typeTable);
		Button btAddNewRow = new Button("New Value");
		btAddNewRow.setDescription("Adds new Value to the table");
		btAddNewRow.setIcon(new ThemeResource("images/edit-add-2.png"));
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
				defaultValue2.addItem(name);
			}
		});
		vLayoutTypetable.addComponent(btAddNewRow);
	}
	
	private void fillTypeTableWithData(final Name[] types) {
		for(int i = 0; i < types.length; i++) {
			typeTable.addItem(new Object[] {types[i].getID(), types[i].getDisplayName(), getDeleteButton(types[i])}, types[i]);
		}
	}
	
	private void addTypeTable() {
		if(this.getComponentArea(vLayoutTypetable) == null) {
			int i = 0;
			while(!this.getComponent(1, i).equals(type)) {
				i++;
			}
			this.insertRow(i+1);
			this.addComponent(vLayoutTypetable, 1, i+1);
		}
	}
	
	private void removeTypeTable() {
		this.removeComponent(vLayoutTypetable);
	}
	
	/**
	 * Creates SADL parameter from tool editor parameter
	 * @return SADL parameter
	 */
	public SADLDescription.Parameter getSadlParameter() {
		SADLDescription.Parameter parameter = new SADLDescription.Parameter(Name.createName(id.getValue(), 
				name.getValue()), (ParameterType) type.getValue(), 
				getEnumList(), getValueOrNull(minValue.getValue()), 
				getValueOrNull(maxValue.getValue()), getValueOrNull(getDefaultValue()), 
				getValueOrNull(description.getValue()));
		parameter.setOptional(optional.getValue());
		return parameter;
	}
	
	private String getDefaultValue() {
		if(type.getValue().equals(ParameterType.ENUM)) {
			return (defaultValue2.getValue() != null ? ((Name) defaultValue2.getValue()).getID() : null);
		} else if(type.getValue().equals(ParameterType.PERCENT)){
			return String.valueOf(defaultSlider.getValue());
		} else {
			return defaultValue.getValue();
		}
	}
	
	private Name[] getEnumList() {
		if (typeTable == null || typeTable.getItemIds().isEmpty())
			return null;
		Name[] names = typeTable.getItemIds().toArray(new Name[0]);
		ArrayList<Name> list = new ArrayList<Name>();
		for (Name name : names) {
			if (!name.getID().isEmpty()) {
				if ("".equals(name.getDisplayName())) {
					name.setDisplayName(null);
				}
				
				list.add(name);
			} else {
				new Notification("Not all ENUM types were generated to text, because id was empty",  Type.WARNING_MESSAGE).show(Page.getCurrent());
			}
		}
		return list.toArray(new Name[0]);
	}
	
	private Button getDeleteButton(final Object itemId) {
		Button btDelete = new Button();
		btDelete.setIcon(new ThemeResource("images/close.png"));
		btDelete.setStyleName(BaseTheme.BUTTON_LINK);
		btDelete.addClickListener(new ClickListener() {
			private static final long serialVersionUID = -3695725710938486562L;

			@Override
			public void buttonClick(ClickEvent event) {
				typeTable.removeItem(itemId);
				defaultValue2.removeItem(itemId);
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
		return type.getValue().toString();
	}
	
	private void removeMinMaxFields() {
		if(this.getComponentArea(maxValue) != null && this.getComponentArea(minValue) != null) {
			this.removeRow(this.getComponentArea(maxValue).getRow1());
			this.removeRow(this.getComponentArea(minValue).getRow1());
		}
	}
	
	private void AddMinMaxFields() {
		if(this.getComponentArea(minValue) == null && this.getComponentArea(maxValue) == null) {
			int i = 0;
			while(!this.getComponent(1, i).equals(type)) {
				i++;
			}
			this.insertRow(i+1);
			this.addComponent(lbMaxValue, 0, i+1);
			this.addComponent(maxValue, 1, i+1);
			this.insertRow(i+2);
			this.addComponent(lbMinValue, 0, i+2);
			this.addComponent(minValue, 1, i+2);
		}
	}
	
	private void changeDefautlToDropDown() {
		if(this.getComponentArea(defaultValue) != null)
			this.replaceComponent(defaultValue, defaultValue2);
		if(this.getComponentArea(defaultSlider) != null)
			this.replaceComponent(defaultSlider, defaultValue2);
	}
	
	private void changeDefaultToTextField()
	{
		if(this.getComponentArea(defaultValue2) != null)
			this.replaceComponent(defaultValue2, defaultValue);
		if(this.getComponentArea(defaultSlider) != null)
			this.replaceComponent(defaultSlider, defaultValue);
	}
	
	private void ChangeDefaultToSlider() {
		if(this.getComponentArea(defaultValue) != null)
			this.replaceComponent(defaultValue, defaultSlider);
		if(this.getComponentArea(defaultValue2) != null)
			this.replaceComponent(defaultValue2, defaultSlider);
	}
	
	private void fillDefaultValueDropdownWithTypeData() {
		@SuppressWarnings("unchecked")
		Collection<Name> ids = (Collection<Name>) typeTable.getItemIds();
		for(Name name : ids) {
			defaultValue2.addItem(name);
			defaultValue2.setItemCaption(name, name.getDisplayName());
		}
	}
	
	@Override
	public String toString() {
		return getValue((id.getValue() != null && !id.getValue().isEmpty() ? id.getValue() :name.getValue())) + " " + getValue(type.getValue().toString());
	}
	
	public boolean validateFields() {
		boolean isValid = maxValue.isValid() && minValue.isValid() && id.isValid();
		return (type.getValue().equals(ParameterType.ENUM) && type.getValue().equals(ParameterType.PERCENT) 
				? isValid :(isValid && defaultValue.isValid()));
	}
}
