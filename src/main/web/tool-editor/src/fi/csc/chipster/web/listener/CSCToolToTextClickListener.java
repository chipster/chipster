package fi.csc.chipster.web.listener;

import com.vaadin.server.Page;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.Notification.Type;

import fi.csc.chipster.web.model.Tool;
import fi.csc.chipster.web.tooledit.ToolEditorUI;
import fi.csc.microarray.description.SADLDescription;

/**
 * Creating text from tool editor
 * @author Gintare Pacauskaite
 *
 */
public class CSCToolToTextClickListener implements ClickListener{

	private static final long serialVersionUID = 1833698097187666721L;
	
	private ToolEditorUI root;
	
	public CSCToolToTextClickListener(ToolEditorUI root) {
		this.root = root;
	}

	@Override
	public void buttonClick(ClickEvent event) {
		Tool tool = root.getToolEditor().getTool();
		if(tool == null) {
			new Notification("Tool is missing, please insert Tool", Type.WARNING_MESSAGE).show(Page.getCurrent());
			return;
		}
		SADLDescription sadlDescription = tool.createSadlDescription();
		if(sadlDescription == null) {
			new Notification("Tool's elements are empty, please fill it up", Type.WARNING_MESSAGE).show(Page.getCurrent());
			return;
		}
		try {
			sadlDescription.addInputs(root.getToolEditor().getSADLInputs());
			sadlDescription.addOutputs(root.getToolEditor().getSADLOutputs());
			root.getToolEditor().addParameters(sadlDescription);
			root.getToolEditor().setHeaderToTextEditor(sadlDescription.toString());
		} catch(Exception e) {
			Notification.show("All elements must be filled up correctly\n" + e.getMessage(), Type.WARNING_MESSAGE);
			e.printStackTrace();
			return;
		}
		
		
	}

}
