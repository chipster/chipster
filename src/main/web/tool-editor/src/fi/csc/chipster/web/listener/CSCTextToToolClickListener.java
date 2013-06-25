package fi.csc.chipster.web.listener;

import java.util.List;

import com.vaadin.server.Page;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;

import fi.csc.chipster.web.tooledit.ToolEditorUI;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class CSCTextToToolClickListener implements ClickListener{
	
	private static final long serialVersionUID = 8517086427914209637L;
	
	private ToolEditorUI root;
	
	public CSCTextToToolClickListener(ToolEditorUI root) {
		this.root = root;
	}

	@Override
	public void buttonClick(ClickEvent event) {
		String text = root.getTextEditor().takeHeader(root.getTextEditor().getText());
		
		ChipsterSADLParser parser = new ChipsterSADLParser();
		try {
			SADLDescription description = parser.parse(text);
			root.getToolEditor().removeItems();
			root.getTreeToolEditor().removeAllChildren();
			root.getToolEditor().addTool(description);
			List<Input> inputs = description.inputs();
			for(int i = 0 ; i < inputs.size(); i++) {
				root.getToolEditor().addInput(inputs.get(i));
			}
			List<Output> outputs = description.outputs();
			for(int i = 0 ; i < outputs.size(); i++) {
				root.getToolEditor().addOutput(outputs.get(i));
			}
			List<Parameter> parameters = description.parameters();
			for(int i = 0 ; i < parameters.size(); i++) {
				root.getToolEditor().addParameter(parameters.get(i));
			}
			
		} catch (Exception e) {
			new Notification("Please fill text area", Type.WARNING_MESSAGE).show(Page.getCurrent());
		}
		
	}

	

}
