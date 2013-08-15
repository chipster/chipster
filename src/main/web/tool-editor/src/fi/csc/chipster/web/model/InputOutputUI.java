package fi.csc.chipster.web.model;

import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.ComboBox;
import com.vaadin.ui.Label;

/**
 * Abstract class for input and output models.
 * @author Gintare Pacauskaite
 *
 */
public abstract class InputOutputUI extends BasicModel {
	
	private static final long serialVersionUID = 7877139214264349071L;
	
	protected Label lbType;
	protected Label lbType2;
	protected Label lbMeta;
	
	protected ComboBox type;
	protected ComboBox type2;
	protected CheckBox cbMeta;
	
	public InputOutputUI() {
	}
	
	protected void initElements() {
		type2 = new ComboBox();
		type2.setImmediate(true);
		type2.setWidth(WIDTH);
		lbMeta = new Label("Meta:");
		cbMeta = new CheckBox();
		cbMeta.setDescription("Is this element Meta data");
		
		optional.setWidth(WIDTH);
		type = new ComboBox();
		type.setWidth(WIDTH);
		type.setNullSelectionAllowed(false);
		type.setImmediate(true);
		type.addItem(SINGLE_FILE);
		type.addItem(MULTI_FILE);
		type.select(SINGLE_FILE);
		type.addValueChangeListener(new ValueChangeListener() {
			private static final long serialVersionUID = -1134955257251483403L;

			@Override
			public void valueChange(ValueChangeEvent event) {
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
	
	@Override
	protected String getType() {
		return type2.getItemCaption(type2.getValue());
	}

}
