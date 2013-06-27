package fi.csc.chipster.web.listener;

import com.vaadin.event.FieldEvents.TextChangeEvent;
import com.vaadin.event.FieldEvents.TextChangeListener;

import fi.csc.chipster.web.model.BasicModel;
import fi.csc.chipster.web.model.Tool;

@SuppressWarnings("serial")
public class CSCTextChangeListener implements TextChangeListener{
	
	private BasicModel model;
	private boolean isId = false;
	
	public CSCTextChangeListener(BasicModel model, boolean isId) {
		this.model = model;
		this.isId = isId;
	}
	
	public CSCTextChangeListener(BasicModel model) {
		this.model = model;
	}

	@Override
	public void textChange(TextChangeEvent event) {
		String text = "";
		if(model instanceof Tool) {
			System.out.println(event.getText());
			if(isId) {
				text = ((Tool) model).getTitleId(event.getText());
			} else {
				text = ((Tool) model).getTitleDisplayName(event.getText());
			}
			model.getToolEditor().getToolEditorUI().getTreeToolEditor().updateToolTitle(text);
		} else {
			text = getValue(event.getText()) + " " + getValue(model.getTypeValue()) + (model.isOptional() ? " OPTIONAL" : "");
			System.out.println("textchange: " + text);
			model.getToolEditor().getToolEditorUI().getTreeToolEditor().setItemCaption(model, text);
		}
		model.setTitleDescriptionValue(text);
		
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
