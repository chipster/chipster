package fi.csc.chipster.web.model;

import com.vaadin.ui.CheckBox;
import com.vaadin.ui.ComboBox;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;

public class Output extends BasicModel{

	private CheckBox optional;
	private ComboBox type;
	
	public GridLayout createOutputUI() {
//		grid.addComponent(new Label("Output"), 0, 0);
		initElements();
		grid.addComponent(type, 0, 1);
		grid.addComponent(optional, 1, 1);
		return grid;
	}
	
	private void initElements() {
		optional = new CheckBox("Optional");
		type = new ComboBox("Type");
		type.setWidth("100px");
		type.addItem("Single file");
	}

}
