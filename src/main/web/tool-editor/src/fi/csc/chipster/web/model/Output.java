package fi.csc.chipster.web.model;

import com.vaadin.ui.ComboBox;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextField;

public class Output extends BasicModel{

	private Label lbOptional;
	private Label lbType;
	private Label lbType2;
	
	private ComboBox optional;
	private ComboBox type;
	private TextField type2;
	
	@Override
	public GridLayout createUI() {
//		grid.addComponent(new Label("Output"), 0, 0);
		initElements();
		addRow(lbType2, type2);
		addRow(lbType, type);
		addRow(lbOptional, optional);
		addRow(lbDescription, description);
		return grid;
	}
	
	private void initElements() {
		lbId.setValue("Output file:");
		lbOptional = new Label("Output is:");
		lbType = new Label("Output is:");
		lbType2 = new Label("Type:");
		
		type2 = new TextField();
		
		optional = new ComboBox();
		optional.addItem(NOT_OPTIONAL);
		optional.addItem(OPTIONAL);
		optional.setNullSelectionAllowed(false);
		type = new ComboBox();
		type.setWidth("100px");
		type.addItem("Single file");
	}
	
	public GridLayout createUIWithData(fi.csc.microarray.description.SADLDescription.Output output) {
		createUI();
		fillWithData(output);
		return grid;
	}
	
	private void fillWithData(fi.csc.microarray.description.SADLDescription.Output output) {
		id.setValue(output.getName().getID());
		name.setValue(output.getName().getDisplayName());
		description.setValue(getValue(output.getComment()));
		if (output.isOptional()) {
			optional.select(OPTIONAL);
		} else {
			optional.select(NOT_OPTIONAL);
		}
	}

}
