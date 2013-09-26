package fi.csc.chipster.web.tooledit;

import com.vaadin.annotations.Theme;
import com.vaadin.server.ThemeResource;
import com.vaadin.server.VaadinRequest;
import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Alignment;
import com.vaadin.ui.Button;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.HorizontalSplitPanel;
import com.vaadin.ui.Label;
import com.vaadin.ui.Panel;
import com.vaadin.ui.UI;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.VerticalSplitPanel;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;

import fi.csc.chipster.web.listener.CSCTextToToolClickListener;
import fi.csc.chipster.web.listener.CSCToolToTextClickListener;

/**
 * Main UI class
 * @author Gintare Pacauskaite
 */
@SuppressWarnings("serial")
@Theme("tool_editor")
public class ToolEditorUI extends UI {
	
	private VerticalSplitPanel vSplitPanel;
	private ToolEditor toolEditor;
	private TextEditor textEditor;
	private TreeToolEditor treeToolEditor;
	private Button btUpdateToolEditor;
	private Button btUpdateTextEditor;


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

        HorizontalLayout buttonPanel = getButtonPanel();
        vvLayout.addComponent(buttonPanel);
        vvLayout.setComponentAlignment(buttonPanel, Alignment.MIDDLE_CENTER);

        vvLayout.setExpandRatio(hSplitpPanel, 5);
        vvLayout.setComponentAlignment(hSplitpPanel, Alignment.TOP_LEFT);
        vvLayout.setMargin(false);
        vvLayout.setSpacing(false);
        hSplitpPanel.setFirstComponent(treeToolEditor);
        hSplitpPanel.setSecondComponent(toolEditor);
        vSplitPanel.setFirstComponent(vvLayout);
        vSplitPanel.setSecondComponent(textEditor);
        hSplitpPanel.setStyleName("topborder");
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
    
    private HorizontalLayout getButtonPanel() {
		HorizontalLayout hLayout = new HorizontalLayout();
		hLayout.setSpacing(true);
		
		btUpdateTextEditor = new Button();
		btUpdateTextEditor.setDescription("Update text area");
		btUpdateTextEditor.setIcon(new ThemeResource("images/arrow_down.png"));
		btUpdateTextEditor.addClickListener(new CSCToolToTextClickListener(this));
		hLayout.addComponent(btUpdateTextEditor);
		btUpdateToolEditor = new Button();
		btUpdateToolEditor.setDescription("Update tool elements");
		btUpdateToolEditor.setIcon(new ThemeResource("images/arrow_up.png"));
		btUpdateToolEditor.addClickListener(new CSCTextToToolClickListener(this));
		hLayout.addComponent(btUpdateToolEditor);
		Button btClearAll = new Button("Clear All");
		btClearAll.addClickListener(new ClickListener() {
			private static final long serialVersionUID = 1487893808578560989L;

			@Override
			public void buttonClick(ClickEvent event) {
				
				ToolEditorUI.this.addWindow(new ConfirmClearAll(ToolEditorUI.this));
			}
		});
		hLayout.addComponent(btClearAll);
		return hLayout;

    }
    
}