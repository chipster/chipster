package fi.csc.chipster.web.tooledit;

import com.vaadin.annotations.Theme;
import com.vaadin.server.VaadinRequest;
import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Alignment;
import com.vaadin.ui.HorizontalSplitPanel;
import com.vaadin.ui.Label;
import com.vaadin.ui.Panel;
import com.vaadin.ui.UI;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.VerticalSplitPanel;

/**
 * Main UI class
 */
@SuppressWarnings("serial")
@Theme("tool_editor")
public class ToolEditorUI extends UI {
	
	private VerticalSplitPanel vSplitPanel;
	private ToolEditor toolEditor;
	private TextEditor textEditor;
	private TreeToolEditor treeToolEditor;

    @Override
    protected void init(VaadinRequest request) {
    	
    	treeToolEditor = new TreeToolEditor(this);
    	toolEditor = new ToolEditor(this);
    	textEditor = new TextEditor(this);
    	final Panel vLayout = new Panel();
    	vSplitPanel = new VerticalSplitPanel();
    	vSplitPanel.setSplitPosition(50, Unit.PERCENTAGE);
    	vSplitPanel.setImmediate(true);
    	vSplitPanel.setLocked(false);
    	vSplitPanel.setWidth("100%");
//    	vSplitPanel.setCaption("Tool Editor2");
    	vLayout.setContent(vSplitPanel);
        setContent(vSplitPanel);
        VerticalLayout vvLayout = new VerticalLayout();
        vvLayout.setSizeFull();
        Label title = new Label("<h2><b>&nbsp;Tool Editor</b></h2>", ContentMode.HTML);
        vvLayout.addComponent(title);
        vvLayout.setComponentAlignment(title, Alignment.TOP_LEFT);
        HorizontalSplitPanel hSplitpPanel = new HorizontalSplitPanel();
        hSplitpPanel.setSizeFull();
        vvLayout.addComponent(hSplitpPanel);
        vvLayout.setExpandRatio(hSplitpPanel, 5);
        vvLayout.setComponentAlignment(hSplitpPanel, Alignment.TOP_LEFT);
        vvLayout.setMargin(false);
        vvLayout.setSpacing(false);
//        vvLayout.setComponentAlignment(hSplitpPanel, Alignment.TOP_LEFT);
        hSplitpPanel.setFirstComponent(treeToolEditor);
        hSplitpPanel.setSecondComponent(toolEditor);
        vSplitPanel.setFirstComponent(vvLayout);
        vSplitPanel.setSecondComponent(textEditor);
    }
    
    public ToolEditor getToolEditor() {
    	return toolEditor;
    }
    
    public TextEditor getTextEditor() {
    	return textEditor;
    }
    
    public TreeToolEditor getTreeToolEditor() {
    	return treeToolEditor;
    }
}