package fi.csc.chipster.web.model;

import java.util.ArrayList;
import java.util.Collection;

import com.vaadin.data.Container;
import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.data.Validator;
import com.vaadin.data.Validator.InvalidValueException;
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

import fi.csc.chipster.web.tooledit.ToolEditor;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLSyntax.ParameterType;

/**
 * Tool editor parameter model
 * 
 * @author Gintare Pacauskaite
 * 
 */
public class Parameter extends BasicModel {

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
	// for enum type
	private ComboBox enumDefaultValue;
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

	@SuppressWarnings("serial")
	private void initElements() {

		lbType = new Label("Type:");
		lbOptional = new Label("Optional:");
		lbDefaultValue = new Label("Default:");
		lbMaxValue = new Label("Maximum value:");
		lbMinValue = new Label("Minimum value:");
		lbId.setValue("Id:");

		// type
		type = new ComboBox();
		type.setImmediate(true);
		type.setWidth(WIDTH);
		for (ParameterType parameterType : ParameterType.values()) {
			type.addItem(parameterType);
		}
		type.setNullSelectionAllowed(false);
		type.select(type.getItemIds().iterator().next());
		type.addValueChangeListener(new ValueChangeListener() {

			@Override
			public void valueChange(ValueChangeEvent event) {
				defaultValue.setValue("");
				minValue.setValue("");
				maxValue.setValue("");
				if (type.getValue().equals(ParameterType.ENUM)) {
					addTypeTable();
					changeDefautlToDropDown();
					fillDefaultValueDropdownWithTypeData();
					removeMinMaxFields();
				} else {
					removeTypeTable();
					changeDefaultToTextField();
					addMinMaxFields();
					if (type.getValue().equals(ParameterType.STRING)
							|| type.getValue().equals(ParameterType.COLUMN_SEL)
							|| type.getValue().equals(ParameterType.INPUT_SEL)
							|| type.getValue().equals(
									ParameterType.METACOLUMN_SEL)
							|| type.getValue().equals(ParameterType.PERCENT)) {
						removeMinMaxFields();
					}
					if (type.getValue().equals(ParameterType.PERCENT)) {
						ChangeDefaultToSlider();
					}
				}
				String text = getValue(getNameValue()) + " "
						+ getValue(getTypeValue());
				setTitleDescriptionValue(text.trim());
				root.getToolEditorUI()
						.getTreeToolEditor()
						.setItemCaption(Parameter.this,
								Parameter.this.toString());
			}
		});

		// default
		defaultValue = new TextField();
		defaultValue.setDescription("The default value of the parameter. Can be empty");
		defaultValue.setImmediate(true);
		defaultValue.addValidator(getDefaultValueValidator());
		defaultValue.setWidth(WIDTH);

		// default for enum
		enumDefaultValue = new ComboBox();
		enumDefaultValue.setImmediate(true);
		enumDefaultValue.setWidth(WIDTH);

		// default for slider
		defaultSlider = new Slider();
		defaultSlider.setWidth(WIDTH);
		defaultSlider.setImmediate(true);
		defaultSlider.setOrientation(SliderOrientation.HORIZONTAL);

		// max value
		maxValue = new TextField();
		maxValue.setDescription("The maximum value of the parameter. Can be empty");
		maxValue.setImmediate(true);
		maxValue.addValidator(getMaxValueValidator());
		maxValue.setWidth(WIDTH);

		// max value
		minValue = new TextField();
		minValue.setWidth(WIDTH);
		minValue.setImmediate(true);
		minValue.setDescription("Minimum value of parameter. Can be empty");
		minValue.addValidator(getMinValueValidator());
		initTypeTable();
	}


	/**
	 * Creats parameter from SADL parameter
	 * 
	 * @param parameter
	 *            SADL parameter
	 * @return tool editor parameter
	 */
	public Parameter createUIWithData(SADLDescription.Parameter parameter) {
		fillWithData(parameter);
		return this;
	}

	/**
	 * Fills up parameter fields with SADL parameter data
	 * 
	 * @param parameter
	 *            SADL parameter
	 */
	public void fillWithData(SADLDescription.Parameter parameter) {
		type.select(parameter.getType());

		if (parameter.getSelectionOptions() != null) {
			fillTypeTableWithData(parameter.getSelectionOptions());
			addTypeTable();
			fillDefaultValueDropdownWithTypeData();

		}
		id.setValue(getValue(parameter.getName().getID()));
		name.setValue(getValue(parameter.getName().getDisplayName()));
		description.setValue(getValue(parameter.getDescription()));
		optional.setValue(parameter.isOptional());
		if (parameter.getDefaultValues().length >= 1) {
			String defaultVal = parameter.getDefaultValues()[0];
			if (parameter.getType().equals(ParameterType.ENUM)) {
				Name itemId = getdefaultObject(defaultVal);
				if (name != null) {
					enumDefaultValue.select(itemId);
				}
			} else if (parameter.getType().equals(ParameterType.PERCENT)) {
				try {
					defaultSlider.setValue(Double.valueOf(defaultVal));
				} catch (NumberFormatException e) {
				}
			} else {
				defaultValue.setValue(defaultVal);
			}
		}
		maxValue.setValue(getValue(parameter.getTo()));
		minValue.setValue(getValue(parameter.getFrom()));
	}

	@SuppressWarnings("unchecked")
	private Name getdefaultObject(String id) {
		for (Name name : (Collection<Name>) enumDefaultValue.getItemIds()) {
			if (name.getID().equals(id)) {
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
			public Field<?> createField(Container container,
					final Object itemId, final Object propertyId,
					Component uiContext) {

				final TextField txtField = new TextField();
				if (propertyId.equals(COLUMN_ID)) {
					txtField.setValue(((Name) itemId).getID());
					typeTable.getContainerProperty(itemId, COLUMN_ID).setValue(
							getValue(((Name) itemId).getID()));
				} else if (propertyId.equals(COLUMN_NAME)) {
					txtField.setValue(((Name) itemId).getDisplayName());
					typeTable.getContainerProperty(itemId, COLUMN_NAME)
							.setValue(
									getValue(((Name) itemId).getDisplayName()));
				}

				txtField.addValueChangeListener(new ValueChangeListener() {
					private static final long serialVersionUID = 7785711696111601557L;

					@Override
					public void valueChange(ValueChangeEvent event) {
						if (propertyId.equals(COLUMN_ID)) {
							((Name) itemId).setID(txtField.getValue());
						} else if (propertyId.equals(COLUMN_NAME)) {
							((Name) itemId).setDisplayName(txtField.getValue());
							enumDefaultValue.setItemCaption(itemId,
									txtField.getValue());
						}
					}
				});

				if (typeTable.getContainerDataSource().size() > 5) {
					typeTable.setPageLength(5);
				} else {
					typeTable.setPageLength(typeTable.getContainerDataSource()
							.size());
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
				DataBoundTransferable transferable = (DataBoundTransferable) event
						.getTransferable();
				Object sourceRow = transferable.getItemId();

				AbstractSelect.AbstractSelectTargetDetails dropData = ((AbstractSelect.AbstractSelectTargetDetails) event
						.getTargetDetails());
				Object targetRow = dropData.getItemIdOver();

				if ((sourceRow == targetRow) || (targetRow == null)) {
					return;
				}
				typeTable.removeItem(sourceRow);
				if (dropData.getDropLocation() == VerticalDropLocation.BOTTOM) {
					typeTable.addItemAfter(targetRow, sourceRow);
				} else {
					Object rowAbove = typeTable.prevItemId(targetRow);
					typeTable.addItemAfter(rowAbove, sourceRow);
				}
				typeTable.getContainerProperty(sourceRow, COLUMN_ACTION)
						.setValue(getDeleteButton(sourceRow));
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
				typeTable.getContainerProperty(name, COLUMN_ACTION).setValue(
						getDeleteButton(name));
				typeTable.setValue(name);
				typeTable.setCurrentPageFirstItemId(name);
				enumDefaultValue.addItem(name);
			}
		});
		vLayoutTypetable.addComponent(btAddNewRow);
	}

	private void fillTypeTableWithData(final Name[] types) {
		for (int i = 0; i < types.length; i++) {
			typeTable.addItem(
					new Object[] { types[i].getID(), types[i].getDisplayName(),
							getDeleteButton(types[i]) }, types[i]);
		}
	}

	private void addTypeTable() {
		if (this.getComponentArea(vLayoutTypetable) == null) {
			int i = 0;
			while (!this.getComponent(1, i).equals(type)) {
				i++;
			}
			this.insertRow(i + 1);
			this.addComponent(vLayoutTypetable, 1, i + 1);
		}
	}

	private void removeTypeTable() {
		this.removeComponent(vLayoutTypetable);
	}

	/**
	 * Creates SADL parameter from tool editor parameter
	 * 
	 * @return SADL parameter
	 */
	public SADLDescription.Parameter getSadlParameter() {
		SADLDescription.Parameter parameter = new SADLDescription.Parameter(
				Name.createName(id.getValue(), name.getValue()),
				(ParameterType) type.getValue(), getEnumList(),
				getValueOrNull(minValue.getValue()),
				getValueOrNull(maxValue.getValue()),
				getValueOrNull(getDefaultValue()),
				getValueOrNull(description.getValue()));
		parameter.setOptional(optional.getValue());
		return parameter;
	}

	private String getDefaultValue() {
		if (type.getValue().equals(ParameterType.ENUM)) {
			return (enumDefaultValue.getValue() != null ? ((Name) enumDefaultValue
					.getValue()).getID() : null);
		} else if (type.getValue().equals(ParameterType.PERCENT)) {
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
				new Notification(
						"Not all ENUM types were generated to text, because id was empty",
						Type.WARNING_MESSAGE).show(Page.getCurrent());
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
				enumDefaultValue.removeItem(itemId);
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
		if (this.getComponentArea(maxValue) != null
				&& this.getComponentArea(minValue) != null) {
			this.removeRow(this.getComponentArea(maxValue).getRow1());
			this.removeRow(this.getComponentArea(minValue).getRow1());
		}
	}

	private void addMinMaxFields() {
		if (this.getComponentArea(minValue) == null
				&& this.getComponentArea(maxValue) == null) {
			int i = 0;
			while (!this.getComponent(1, i).equals(type)) {
				i++;
			}
			this.insertRow(i + 1);
			this.addComponent(lbMaxValue, 0, i + 1);
			this.addComponent(maxValue, 1, i + 1);
			this.insertRow(i + 2);
			this.addComponent(lbMinValue, 0, i + 2);
			this.addComponent(minValue, 1, i + 2);
		}
	}

	private void changeDefautlToDropDown() {
		if (this.getComponentArea(defaultValue) != null)
			this.replaceComponent(defaultValue, enumDefaultValue);
		if (this.getComponentArea(defaultSlider) != null)
			this.replaceComponent(defaultSlider, enumDefaultValue);
	}

	private void changeDefaultToTextField() {
		if (this.getComponentArea(enumDefaultValue) != null)
			this.replaceComponent(enumDefaultValue, defaultValue);
		if (this.getComponentArea(defaultSlider) != null)
			this.replaceComponent(defaultSlider, defaultValue);
	}

	private void ChangeDefaultToSlider() {
		if (this.getComponentArea(defaultValue) != null)
			this.replaceComponent(defaultValue, defaultSlider);
		if (this.getComponentArea(enumDefaultValue) != null)
			this.replaceComponent(enumDefaultValue, defaultSlider);
	}

	private void fillDefaultValueDropdownWithTypeData() {
		@SuppressWarnings("unchecked")
		Collection<Name> ids = (Collection<Name>) typeTable.getItemIds();
		for (Name name : ids) {
			enumDefaultValue.addItem(name);
			enumDefaultValue.setItemCaption(name, name.getDisplayName());
		}
	}

	@Override
	public String toString() {
		return getValue((id.getValue() != null && !id.getValue().isEmpty() ? id
				.getValue() : name.getValue()))
				+ " "
				+ getValue(type.getValue().toString());
	}

	public boolean validateFields() {
		boolean isValid = maxValue.isValid() && minValue.isValid()
				&& id.isValid();
		return (type.getValue().equals(ParameterType.ENUM)
				&& type.getValue().equals(ParameterType.PERCENT) ? isValid
				: (isValid && defaultValue.isValid()));
	}


	

	@SuppressWarnings("serial")
	private Validator getDefaultValueValidator() {
		return new Validator() {

			@Override
			public void validate(Object value) throws InvalidValueException {
				
				String v = (String) value;
				
				// empty is ok
				if (v.isEmpty()) {
					return;
				}

				ParameterType t = (ParameterType) type.getValue();

				if (t == ParameterType.INTEGER) {
					// must be integer
					int defaultValueInt = getInt(v);
					
					// if min value exists, must be equal or greater
					isGreaterOrEqualToMinValueInt(defaultValueInt);
					
					// if max value exist, must be smaller or equal
					isSmallerOrEqualToMaxValueInt(defaultValueInt);
				}

				// DECIMAL
				else if (t == ParameterType.DECIMAL) {
					// must be double
					double defaultValueDouble = getDouble(v);
					
					// if min value exists, must be equal or greater
					isGreaterOrEqualToMinValueDouble(defaultValueDouble);

					// if max value exist, must be smaller or equal
					isSmallerOrEqualToMaxValueDouble(defaultValueDouble);
				}
			}
		};
	}

	@SuppressWarnings("serial")
	private Validator getMaxValueValidator() {
		return new Validator() {

			@Override
			public void validate(Object value) throws InvalidValueException {
				String v = (String) value;

				// empty is ok
				if (v.isEmpty()) {
					return;
				}

				// integer
				ParameterType t = (ParameterType) type.getValue();
				if (t == ParameterType.INTEGER) {
					// must be integer
					int maxValueInt = getInt(v);

					// if min value exists, must be equal or greater
					isGreaterOrEqualToMinValueInt(maxValueInt);

					// TODO check default after max change
				} 

				
				// DECIMAL
				else if (t == ParameterType.DECIMAL) {
					double maxValueDouble = getDouble(v);
					isGreaterOrEqualToMinValueDouble(maxValueDouble);
				}
			}
		};
	}

	
	@SuppressWarnings("serial")
	private Validator getMinValueValidator() {
		return new Validator() {

			@Override
			public void validate(Object value) throws InvalidValueException {
				String v = (String) value;

				// empty is ok
				if (v.isEmpty()) {
					return;
				}

				// integer
				ParameterType t = (ParameterType) type.getValue();
				if (t == ParameterType.INTEGER) {
					// must be integer
					int minValueInt = getInt(v);

					// if min value exists, must be equal or greater
					isSmallerOrEqualToMaxValueInt(minValueInt);
				} 

				
				// DECIMAL
				else if (t == ParameterType.DECIMAL) {
					double minValueDouble = getDouble(v);
					isSmallerOrEqualToMaxValueDouble(minValueDouble);
				}
			}
		};
	}

	
	/**
	 * @throws InvalidValueException if min value exists and i is smaller than min value
	 * @param i
	 */
	private void isGreaterOrEqualToMinValueInt(int i) {
		if (!minValue.getValue().isEmpty()) {
			try {
				int minValueInt = Integer.parseInt(minValue.getValue());
				if (i < minValueInt) {
					throw new InvalidValueException("Smaller than minimum value");
				}
			} catch (NumberFormatException e) {
				// ignore rubbish min value
			}
		}
	}

	/**
	 * @throws InvalidValueException if max value exists and i is greater than max value
	 * @param i
	 */
	private void isSmallerOrEqualToMaxValueInt(int i) {
		if (!maxValue.getValue().isEmpty()) {
			try {
				int maxValueInt = Integer.parseInt(maxValue.getValue());
				if (i > maxValueInt) {
					throw new InvalidValueException("Greater than maximum value");
				}
			} catch (NumberFormatException e) {
				// ignore rubbish max value
			}
		}
	}

	/**
	 * @throws InvalidValueException if min value exists and i is smaller than min value
	 * @param i
	 */
	private void isGreaterOrEqualToMinValueDouble(double d) {
		if (!minValue.getValue().isEmpty()) {
			try {
				double minValueDouble = Double.parseDouble(minValue.getValue());
				if (d < minValueDouble) {
					throw new InvalidValueException("Smaller than minimum value");
				}
			} catch (NumberFormatException e) {
				// ignore rubbish min value
			}
		}
	}

	/**
	 * @throws InvalidValueException if max value exists and i is greater than max value
	 * @param i
	 */
	private void isSmallerOrEqualToMaxValueDouble(double d) {
		if (!maxValue.getValue().isEmpty()) {
			try {
				double maxValueDouble = Double.parseDouble(maxValue.getValue());
				if (d > maxValueDouble) {
					throw new InvalidValueException("Greater than maximum value");
				}
			} catch (NumberFormatException e) {
				// ignore rubbish max value
			}
		}
	}

	/**
	 * 
	 * @param s
	 * @return
	 * @throws InvalidValueException if s can not be parsed to int
	 */
	private int getInt(String s) {
		int i;
		try {
			i = Integer.parseInt(s);
		} catch (NumberFormatException e) {
			throw new InvalidValueException("Not an integer");
		}

		return i; 
	}

	
	/**
	 * 
	 * @param s
	 * @return
	 * @throws InvalidValueException if s can not be parsed to double
	 */
	private double getDouble(String s) {
		double d;
		try {
			d = Double.parseDouble(s);
		} catch (NumberFormatException e) {
			throw new InvalidValueException("Not a decimal number");
		}

		return d; 
	}


}
