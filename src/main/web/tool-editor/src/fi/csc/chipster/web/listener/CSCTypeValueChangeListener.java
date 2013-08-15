package fi.csc.chipster.web.listener;

import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;

import fi.csc.chipster.web.model.BasicModel;

/**
 * Changing summary tree node captions and on the right side of tool editor, when type has changed
 * @author Gintare Pacauskaite
 *
 */
public class CSCTypeValueChangeListener implements ValueChangeListener{
	private static final long serialVersionUID = 8903198740449805236L;
	private BasicModel model;
	
	public CSCTypeValueChangeListener(BasicModel model) {
		this.model = model;
	}

	@Override
	public void valueChange(ValueChangeEvent event) {
		String text = getValue(model.getNameValue()) + " " + getValue(model.getTypeValue());
		model.setTitleDescriptionValue(text);
		model.getToolEditor().getToolEditorUI().getTreeToolEditor().setItemCaption(model, text);
		
	}
	
	private String getValue(String text) {
		return (text == null || text.isEmpty() || text.trim().isEmpty() ? "" : text);
	}

}
