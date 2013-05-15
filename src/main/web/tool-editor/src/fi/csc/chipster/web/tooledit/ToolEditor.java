package fi.csc.chipster.web.tooledit;

import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextField;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Button.ClickEvent;

import fi.csc.chipster.web.model.Input;
import fi.csc.chipster.web.model.Output;
import fi.csc.chipster.web.model.Parameter;
import fi.csc.chipster.web.model.Tool;



public class ToolEditor extends VerticalLayout{
	
	public ToolEditor() {
		init();
		this.addComponent(setCaption("Tool", false));
		this.addComponent(new Tool().createToolUI());
		this.addComponent(setCaption("Input", false));
		this.addComponent(new Input().createInputUI());
		this.addComponent(setCaption("Output", false));
		this.addComponent(new Output().createOutputUI());
		this.addComponent(setCaption("Parameter", true));
		this.addComponent(new Parameter().createParameterUI());
	}
	
	private void init() {
		this.setImmediate(true);
	}
	
	private HorizontalLayout setCaption(String name, boolean button) {
		HorizontalLayout h = new HorizontalLayout();
		h.addComponent(new Label("<b>" + name + "</b>", ContentMode.HTML));
		if(button) {
			Button b = new Button("+");
			System.out.println("bee");
			h.addComponent(b);
			h.setSpacing(true);
			b.addClickListener(new ClickListener() {
				
				@Override
				public void buttonClick(ClickEvent event) {
					ToolEditor.this.addComponent(new Parameter().createParameterUI());
					
				}
			});
		}
		return h;
	}

}
