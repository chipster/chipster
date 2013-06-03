package fi.csc.chipster.web.model;

import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.ComboBox;
import com.vaadin.ui.Label;

public abstract class InputOutputUI extends BasicModel {

	protected Label lbOptional;
	protected Label lbType;
	protected Label lbType2;
	protected Label lbMeta;
	
	protected ComboBox optional;
	protected ComboBox type;
	protected ComboBox type2;
	protected CheckBox cbMeta;
	
	public InputOutputUI() {
		// TODO Auto-generated constructor stub
	}
	
protected void initElements() {
		type2 = new ComboBox();
		type2.setImmediate(true);
		type2.setWidth(WIDTH);
		lbMeta = new Label("Meta:");
		cbMeta = new CheckBox();
		
		optional = new ComboBox();
		optional.setWidth(WIDTH);
		optional.setNullSelectionAllowed(false);
		optional.addItem(NOT_OPTIONAL);
		optional.addItem(OPTIONAL);
		optional.select(NOT_OPTIONAL);
		type = new ComboBox();
		type.setWidth(WIDTH);
		type.setNullSelectionAllowed(false);
		type.setImmediate(true);
		type.addItem(SINGLE_FILE);
		type.addItem(MULTI_FILE);
		type.select(SINGLE_FILE);
		type.addValueChangeListener(new ValueChangeListener() {
			@Override
			public void valueChange(ValueChangeEvent event) {
				// TODO Auto-generated method stub
				if(type.getValue().toString().contentEquals(SINGLE_FILE)) {
					getSingleFileUI();
				} else if(type.getValue().toString().contentEquals(MULTI_FILE)){
					getMultipleFilesUI();
				}
			}
		});
	}

	protected void createBasicUI() {
		addRow(lbType2, type2);
		addRow(lbType, type);
		addRow(lbMeta, cbMeta);
		addRow(lbOptional, optional);
		addRow(lbDescription, description);
	}

}
