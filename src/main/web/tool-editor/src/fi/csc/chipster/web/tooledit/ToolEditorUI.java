package fi.csc.chipster.web.tooledit;

import com.vaadin.annotations.Theme;
import com.vaadin.server.Sizeable;
import com.vaadin.server.VaadinRequest;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Panel;
import com.vaadin.ui.UI;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.VerticalSplitPanel;

import fi.csc.chipster.web.model.Tool;

/**
 * Main UI class
 */
@SuppressWarnings("serial")
//@Theme("tool_editor")
public class ToolEditorUI extends UI {
	
	private VerticalSplitPanel vSplitPanel;

    @Override
    protected void init(VaadinRequest request) {
    	
//    	VerticalLayout vertical = new VerticalLayout();
    	
    	final Panel vLayout = new Panel();
    	vSplitPanel = new VerticalSplitPanel();
    	vSplitPanel.setSplitPosition(50, Unit.PERCENTAGE);
    	vSplitPanel.setImmediate(true);
    	vSplitPanel.setLocked(false);
//    	vSplitPanel.setHeight("50%");
    	vSplitPanel.setWidth("100%");
    	vLayout.setContent(vSplitPanel);
        setContent(vSplitPanel);
        
        vSplitPanel.setFirstComponent(new ToolEditor());
        vSplitPanel.setSecondComponent(new TextEditor());
        
//        Button button = new Button("Click Me ");
//        button.addClickListener(new Button.ClickListener() {
//            public void buttonClick(ClickEvent event) {
//                vSplitPanel.addComponent(new Label("Thank you for clicking"));
//            }
//        });
//        vSplitPanel.addComponent(button);
    }
}