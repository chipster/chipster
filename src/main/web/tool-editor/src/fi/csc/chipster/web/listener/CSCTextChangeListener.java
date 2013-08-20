package fi.csc.chipster.web.listener;

import com.vaadin.event.FieldEvents.TextChangeEvent;
import com.vaadin.event.FieldEvents.TextChangeListener;

import fi.csc.chipster.web.model.BasicModel;
import fi.csc.chipster.web.model.Tool;
/**
 * Listener for changing text in summer tree node captions and on the right side of tool editor
 * @author Gintare Pacauskaite
 *
 */
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
			if(isId) {
				text = ((Tool) model).getTitleId(event.getText());
			} else {
				text = ((Tool) model).getTitleDisplayName(event.getText());
			}
			model.getToolEditor().getToolEditorUI().getTreeToolEditor().updateToolTitle(text);
		} else {
			if(!isId && (model.getId() != null && !model.getId().isEmpty())) {
				return;
			}
			if(event.getText().isEmpty() && isId) {
				text = getValue(model.getNameValue()) + " " + getValue(model.getTypeValue());
			} else {
				text = getValue(event.getText()) + " " + getValue(model.getTypeValue());
			}
			model.getToolEditor().getToolEditorUI().getTreeToolEditor().setItemCaption(model, text);
		}
		model.setTitleDescriptionValue(text);
	}
	
	private String getValue(String text) {
		return (text == null || text.isEmpty() || text.trim().isEmpty() ? "" : text);
	}
}
