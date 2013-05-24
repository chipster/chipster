package fi.csc.chipster.web.tooledit;

import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.model.Input;
import fi.csc.chipster.web.model.Output;
import fi.csc.chipster.web.model.Parameter;
import fi.csc.chipster.web.model.Tool;
import fi.csc.microarray.description.SADLDescription;



public class ToolEditor extends VerticalLayout{
	
	public ToolEditor() {
		init();
//		this.addComponent(setCaption("Tool", false));
////		this.addComponent(new Tool().createToolUI());
//		this.addComponent(setCaption("Input", false));
////		this.addComponent(new Input().createInputUI());
//		this.addComponent(setCaption("Output", false));
////		this.addComponent(new Output().createOutputUI());
//		this.addComponent(setCaption("Parameter", true));
//		this.addComponent(new Parameter().createParameterUI());
	}
	
	private void init() {
		this.setImmediate(true);
	}
	
	private HorizontalLayout setCaption(String name, boolean button) {
		HorizontalLayout h = new HorizontalLayout();
		h.addComponent(new Label("<b>" + name + "</b>", ContentMode.HTML));
		if(button) {
			Button b = new Button("+");
			h.addComponent(b);
			h.setSpacing(true);
			b.addClickListener(new ClickListener() {
				
				@Override
				public void buttonClick(ClickEvent event) {
					ToolEditor.this.addComponent(new Parameter().createUI());
					
				}
			});
		}
		return h;
	}
	public void addTool() {
		this.addComponent(new Tool().createUI());
	}
	
	public void addTool(SADLDescription sadlDescription) {
		this.addComponent(new Tool().createUIwithData(sadlDescription));
	}
	
	public void addInput() {
		this.addComponent(new Input().createUI());
	}
	
	public void addInput(fi.csc.microarray.description.SADLDescription.Input input) {
		this.addComponent(new Input().createUIWithData(input));
	}
	
	public void addOutput() {
		this.addComponent(new Output().createUI());
	}
	
	public void addOutput(fi.csc.microarray.description.SADLDescription.Output output) {
		this.addComponent(new Output().createUIWithData(output));
	}
	
	public void addParameter() {
		this.addComponent(new Parameter().createUI());
	}
	
	public void addParameter(fi.csc.microarray.description.SADLDescription.Parameter parameter) {
		this.addComponent(new Parameter().createUIWithData(parameter));
	}


}
