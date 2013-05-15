package fi.csc.chipster.web.model;

import com.vaadin.ui.CheckBox;
import com.vaadin.ui.ComboBox;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextField;

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
	
	public GridLayout createParameterUI() {
		
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
		defaultValue = new TextField();
		optional = new ComboBox();
		optional.addItem("not optional");
		optional.addItem("optional");
		optional.setNullSelectionAllowed(false);
		maxValue = new TextField();
		minValue = new TextField();
	}

}
