package fi.csc.chipster.web.listener;

import com.vaadin.event.FieldEvents.TextChangeEvent;
import com.vaadin.event.FieldEvents.TextChangeListener;

import fi.csc.chipster.web.model.BasicModel;

@SuppressWarnings("serial")
public class CSCTextChangeListener implements TextChangeListener{
	
	private BasicModel model;
	
	public CSCTextChangeListener(BasicModel model) {
		this.model = model;
	}

	@Override
	public void textChange(TextChangeEvent event) {
		String text = getValue(event.getText()) + " " + getValue(model.getTypeValue());
		model.setTitleDescriptionValue(text.trim());
		model.getToolEditor().getToolEditorUI().getTreeToolEditor().setItemCaption(model, text + (model.isOptional() ? " OPTIONAL" : ""));
		
	}
	
	private String getValue(String text) {
		return (text == null || text.isEmpty() || text.trim().isEmpty() ? "" : text);
	}

//	@Override
//	public void valueChange(ValueChangeEvent event) {
//		// TODO Auto-generated method stub
//		String text = getValue(name) + " " + getValue(type);
//		model.setTitleDescriptionValue(text.trim());
//	}

}
