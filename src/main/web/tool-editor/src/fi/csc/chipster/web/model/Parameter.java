package fi.csc.chipster.web.model;

import com.vaadin.ui.CheckBox;
import com.vaadin.ui.ComboBox;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextField;

public class Parameter extends BasicModel{
	
	private ComboBox type;
	private TextField defaultValue;
	private CheckBox optional;
	
	public GridLayout createParameterUI() {
		
//		grid.addComponent(new Label("Parameter"), 0, 0);
		initElements();
		grid.addComponent(type, 0, 1);
		grid.addComponent(defaultValue, 1, 1);
		grid.addComponent(optional, 0, 2);
		return grid;
	}
	
	private void initElements() {
		type = new ComboBox("Type");
		type.setWidth("100px");
		defaultValue = new TextField("Default");
		optional = new CheckBox("Optional");
	}

}
