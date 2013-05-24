package fi.csc.chipster.web.model;

import com.vaadin.ui.ComboBox;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;

import fi.csc.chipster.web.list.Type;

public class Input extends BasicModel{
	
	private Label lbOptional;
	private Label lbType;
	private Label lbType2;
	
	private ComboBox optional;
	private ComboBox type;
	private ComboBox type2;
	
	@Override
	public GridLayout createUI() {
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
		
		type2 = new ComboBox();
		for(Type type : Type.values()) {
			type2.addItem(type.getName());
		}
		
		optional = new ComboBox();
		optional.setNullSelectionAllowed(false);
		optional.addItem(NOT_OPTIONAL);
		optional.addItem(OPTIONAL);
		type = new ComboBox();
		type.setWidth("100px");
		type.addItem("Single file");
	}
	
	public GridLayout createUIWithData(fi.csc.microarray.description.SADLDescription.Input input) {
		createUI();
		fillWithData(input);
		return grid;
	}
	
	private void fillWithData(fi.csc.microarray.description.SADLDescription.Input input) {
		type2.select(input.getType().getName());
		id.setValue(input.getName().getID());
		name.setValue(input.getName().getDisplayName());
		description.setValue(getValue(input.getComment()));
		if (input.isOptional()) {
			optional.select(OPTIONAL);
		} else {
			optional.select(NOT_OPTIONAL);
		}
	}

}
