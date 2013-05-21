package fi.csc.chipster.web.model;

import com.vaadin.ui.ComboBox;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Table;
import com.vaadin.ui.TextField;

import fi.csc.microarray.description.SADLDescription.Name;
import fi.csc.microarray.description.SADLSyntax.ParameterType;


public class Parameter extends BasicModel{
	
	private Label lbType;
	private Label lbOptional;
	private Label lbDefaultValue;
	private Label lbMaxValue;
	private Label lbMinValue;
	
	private ComboBox type;
	private TextField defaultValue;
	private ComboBox optional;
	private TextField maxValue;
	private TextField minValue;
	private Table typeTable;
	
	@Override
	public GridLayout createUI() {
		
//		grid.addComponent(new Label("Parameter"), 0, 0);
		initElements();
		addRow(lbType, type);
		addRow(lbMinValue, minValue);
		addRow(lbDefaultValue, defaultValue);
		addRow(lbMaxValue, maxValue);
		addRow(lbOptional, optional);
		addRow(lbDescription, description);
		return grid;
	}
	
	private void initElements() {
		
		lbType = new Label("Type:");
		lbOptional = new Label("Parameter is:");
		lbDefaultValue = new Label("Default:");
		lbMaxValue = new Label("Maximum value:");
		lbMinValue = new Label("Minimum value:");
		lbId.setValue("User defined parameter");
		
		type = new ComboBox();
		type.setWidth("100px");
		for(ParameterType parameterType : ParameterType.values()) {
			type.addItem(parameterType);
		}
		defaultValue = new TextField();
		optional = new ComboBox();
		optional.addItem(NOT_OPTIONAL);
		optional.addItem(OPTIONAL);
		optional.setNullSelectionAllowed(false);
		maxValue = new TextField();
		minValue = new TextField();
	}
	
	public GridLayout createUIWithData(fi.csc.microarray.description.SADLDescription.Parameter parameter) {
		createUI();
		fillWithData(parameter);
		return grid;
	}
	
	private void fillWithData(fi.csc.microarray.description.SADLDescription.Parameter parameter) {
		parameter.getType();
		type.select(parameter.getType());
		
		if(parameter.getSelectionOptions() != null) {
			initTypeTable();
			fillTypeTableWithData(parameter.getSelectionOptions());
			addTypeTable();
		}
		id.setValue(parameter.getName().getID());
		name.setValue(parameter.getName().getDisplayName());
		description.setValue(parameter.getComment());
		if (parameter.isOptional()) {
			optional.select(OPTIONAL);
		} else {
			optional.select(NOT_OPTIONAL);
		}
		if(parameter.getDefaultValues().length == 1)
			defaultValue.setValue(parameter.getDefaultValue());
	}
	
	private void initTypeTable() {
		typeTable = new Table();
		typeTable.addContainerProperty("ID", String.class, null);
		typeTable.addContainerProperty("Name", String.class, null);
		typeTable.setHeight("100px");
	}
	
	private void fillTypeTableWithData(Name[] types) {
		for(int i = 0; i < types.length; i++) {
			typeTable.addItem(new Object[] {types[i].getID(), types[i].getDisplayName()}, new Integer(i));
			
		}
	}
	
	private void addTypeTable() {
		int i = 0;
		while(!grid.getComponent(1, i).equals(type)) {
			i++;
		}
		grid.insertRow(i+1);
		grid.addComponent(typeTable, 1, i+1);
	}

}
