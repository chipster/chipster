package fi.csc.chipster.web.tooledit;

import java.util.ArrayList;
import java.util.List;

import com.vaadin.server.Page;
import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.Component;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.Notification;
import com.vaadin.ui.Notification.Type;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.model.BasicModel;
import fi.csc.chipster.web.model.Input;
import fi.csc.chipster.web.model.Output;
import fi.csc.chipster.web.model.Parameter;
import fi.csc.chipster.web.model.Tool;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLGenerator;



public class ToolEditor extends VerticalLayout{
	
	private final int INPUT = 0;
	private final int OUTPUT = 1;
	private HorizontalLayout hlHeader;
	private ToolEditorUI root;
	
	public ToolEditor(ToolEditorUI root) {
		this.root = root;
		init();
		buttonHeader();
//		this.addComponent(setCaption("Tool", false));
////		this.addComponent(new Tool().createToolUI());
//		this.addComponent(setCaption("Input", false));
////		this.addComponent(new Input().createInputUI());
//		this.addComponent(setCaption("Output", false));
////		this.addComponent(new Output().createOutputUI());
//		this.addComponent(setCaption("Parameter", true));
//		this.addComponent(new Parameter().createParameterUI());
	}
	
	private void buttonHeader() {
		hlHeader = new HorizontalLayout();
		hlHeader.setImmediate(true);
		Button btAddTool = new Button("Add Tool");
		Button btAddInput = new Button("Add Input");
		Button btAddOutput = new Button("Add Output");
		Button btAddParameter = new Button("Add Parameter");
		Button btUpdate = new Button("Update");
		hlHeader.setSpacing(true);
		hlHeader.setMargin(true);
		hlHeader.addComponent(btAddTool);
		hlHeader.addComponent(btAddInput);
		hlHeader.addComponent(btAddOutput);
		hlHeader.addComponent(btAddParameter);
		hlHeader.addComponent(btUpdate);
		this.addComponent(hlHeader);
		btAddTool.addClickListener(new ClickListener() {
			
			@Override
			public void buttonClick(ClickEvent event) {
				// TODO Auto-generated method stub
				ToolEditor.this.addTool();
			}
		});
		btAddInput.addClickListener(new ClickListener() {
			
			@Override
			public void buttonClick(ClickEvent event) {
				// TODO Auto-generated method stub
				ToolEditor.this.addInput();
			}
		});
		
		btAddOutput.addClickListener(new ClickListener() {
			
			@Override
			public void buttonClick(ClickEvent event) {
				// TODO Auto-generated method stub
				ToolEditor.this.addOutput();
			}
		});
		
		btAddParameter.addClickListener(new ClickListener() {
			
			@Override
			public void buttonClick(ClickEvent event) {
				// TODO Auto-generated method stub
				ToolEditor.this.addParameter();
			}
		});
		
		btUpdate.addClickListener(new ClickListener() {
			
			@Override
			public void buttonClick(ClickEvent event) {
				Tool tool = getTool();
				if(tool == null) {
					new Notification("Tool is missing, please insert Tool", Type.WARNING_MESSAGE).show(Page.getCurrent());
					return;
				}
				SADLDescription sadlDescription = tool.createSadlDescription();
				if(sadlDescription == null) {
					new Notification("Tool's elements are empty, please fill it up", Type.WARNING_MESSAGE).show(Page.getCurrent());
					return;
				}
				sadlDescription.addInputs(getDaslInputs());
				sadlDescription.addOutputs(getDaslOutputs());
				addParameters(sadlDescription);
				System.out.println(sadlDescription.toString());
				setHeaderToTextEditor(sadlDescription.toString());
			}
		});
	}
	
	private void init() {
		this.setImmediate(true);
	}
	
//	private HorizontalLayout setCaption(String name, boolean button) {
//		HorizontalLayout h = new HorizontalLayout();
//		h.addComponent(new Label("<b>" + name + "</b>", ContentMode.HTML));
//		if(button) {
//			Button b = new Button("+");
//			h.addComponent(b);
//			h.setSpacing(true);
//			b.addClickListener(new ClickListener() {
//				
//				@Override
//				public void buttonClick(ClickEvent event) {
//					ToolEditor.this.addComponent(new Parameter().createUI());
//					
//				}
//			});
//		}
//		return h;
//	}
	public void addTool() {
		if(getTool() == null)
			this.addComponent(new Tool(this).createUI(), 1);
	}
	
	public void addTool(SADLDescription sadlDescription) {
		this.addComponent(new Tool(this).createUIwithData(sadlDescription));
	}
	
	public void addInput() {
		this.addComponent(new Input(this).createUI(), getLastIndex(INPUT)+1);
	}
	
	public void addInput(SADLDescription.Input input) {
		this.addComponent(new Input(this).createUIWithData(input));
	}
	
	public void addOutput() {
		this.addComponent(new Output(this).createUI(), getLastIndex(OUTPUT)+1);
	}
	
	public void addOutput(SADLDescription.Output output) {
		this.addComponent(new Output(this).createUIWithData(output));
	}
	
	public void addParameter() {
		this.addComponent(new Parameter(this).createUI());
	}
	
	public void addParameter(SADLDescription.Parameter parameter) {
		this.addComponent(new Parameter(this).createUIWithData(parameter));
	}
	
	public void removeItems() {
		this.removeAllComponents();
		this.addComponent(hlHeader);
	}

	public Tool getTool() {
		for(int i = 0; i < this.getComponentCount(); i++) {
			if(this.getComponent(i) instanceof Tool) {
				return (Tool) this.getComponent(i);
			}
		}
		return null;
	}
	
	public List<SADLDescription.Input> getDaslInputs() {
		ArrayList<SADLDescription.Input> inputs = new ArrayList<SADLDescription.Input>();
		for(int i = 0; i < this.getComponentCount(); i++) {
			if(this.getComponent(i) instanceof Input) {
				inputs.add(((Input) this.getComponent(i)).getSadlInput());
			} else {
				if(!inputs.isEmpty()) {
					break;
				}
			}
		}
		return inputs;
	}
	
	public List<SADLDescription.Output> getDaslOutputs() {
		ArrayList<SADLDescription.Output> outputs = new ArrayList<SADLDescription.Output>();
		for(int i = 0; i < this.getComponentCount(); i++) {
			if(this.getComponent(i) instanceof Output) {
				outputs.add(((Output) this.getComponent(i)).getSadlOutput());
			} else {
				if(!outputs.isEmpty()) {
					break;
				}
			}
		}
		return outputs;
	}
	
	public void addParameters(SADLDescription sadlDescription) {
		boolean wasParameter = false;
		for(int i = 0; i < this.getComponentCount(); i++) {
			if(this.getComponent(i) instanceof Parameter) {
				sadlDescription.addParameter(((Parameter) this.getComponent(i)).getSadlParameter());
				wasParameter = true;
			} else {
				if(wasParameter) {
					break;
				}
			}
		}
//		return sadlDescription;
	}
	
	private void setHeaderToTextEditor(String text) {
		System.out.println("iraso");
		root.getTextEditor().setText(createHeader(text));
	}
	
	private String createHeader(String text) {
		StringBuilder str = new StringBuilder();
		String[] lines = text.split(TextEditor.NEW_LINE);
		for(String line : lines) {
			str.append("# ");
			str.append(line);
			str.append(TextEditor.NEW_LINE);
		}
		System.out.println("grazina");
		return str.toString();
	}
	
	private int getLastIndex(int type) {
		int toolIndex = -1;
		int outputIndex = -1;
		int inputIndex = -1;
		for(int i = 0; i < this.getComponentCount(); i++) {
			if(this.getComponent(i) instanceof Input) {
				inputIndex = i;
			} else if(this.getComponent(i) instanceof Output) {
				outputIndex = i;
			} else if(this.getComponent(i) instanceof Tool) {
				toolIndex = i;
			} else if(this.getComponent(i) instanceof Parameter) {
				break;
			}
		}
		if(type == INPUT) {
			if(inputIndex > -1)
				return inputIndex;
			else if(toolIndex > -1)
				return toolIndex;
			else
				return 0;
		} else if(type == OUTPUT) {
			if(outputIndex > -1)
				return outputIndex;
			else if(inputIndex > -1)
				return inputIndex;
			else if(toolIndex > -1)
				return toolIndex;
			else 
				return 0;
		} else 
			return 0;
	}
	
	public void moveUpComponent(BasicModel component) {
		int index = this.getComponentIndex(component);
		if(index-1 > 0) {
			Component com = this.getComponent(index-1);
			if((component instanceof Input && com instanceof Input) || 
					(component instanceof Output && com instanceof Output) || 
					(component instanceof Parameter && com instanceof Parameter)) {
				this.removeComponent(component);
				this.addComponent(component, index-1);
			}
		}
	}
	
	public void moveDownComponent(BasicModel component) {
		int index = this.getComponentIndex(component);
		if(index+1 < this.getComponentCount()+1) {
			Component com = this.getComponent(index+1);
			if((component instanceof Input && com instanceof Input) || 
					(component instanceof Output && com instanceof Output) || 
					(component instanceof Parameter && com instanceof Parameter)) {
				this.removeComponent(component);
				this.addComponent(component, index+1);
			}
		}
	}
}
