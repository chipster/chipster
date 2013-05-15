package fi.csc.chipster.web.model;

import com.vaadin.ui.ComboBox;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextField;

public class Input extends BasicModel{
	
	private Label lbOptional;
	private Label lbType;
	private Label lbType2;
	
	private ComboBox optional;
	private ComboBox type;
	private TextField type2;
	
	public GridLayout createInputUI() {
//		grid.addComponent(new Label("Input"), 0, 0);
		initElements();
		addRow(lbType2, type2);
		addRow(lbType, type);
		addRow(lbOptional, optional);
		addRow(lbDescription, description);
		
		return grid;
	}
	
	private void initElements() {
		
		lbId.setValue("Input file:");
		lbOptional = new Label("Input is:");
		lbType = new Label("Input is:");
		lbType2 = new Label("Type:");
		
		type2 = new TextField();
		
		optional = new ComboBox();
		optional.addItem("not optional");
		optional.addItem("optional");
		optional.setNullSelectionAllowed(false);
		type = new ComboBox();
		type.setWidth("100px");
		type.addItem("Single file");
	}

}
