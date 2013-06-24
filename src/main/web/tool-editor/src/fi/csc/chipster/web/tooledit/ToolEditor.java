package fi.csc.chipster.web.tooledit;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.vaadin.ui.Component;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.model.BasicModel;
import fi.csc.chipster.web.model.Input;
import fi.csc.chipster.web.model.Output;
import fi.csc.chipster.web.model.Parameter;
import fi.csc.chipster.web.model.Tool;
import fi.csc.microarray.description.SADLDescription;



public class ToolEditor extends VerticalLayout{
	
	private static final long serialVersionUID = -4777411581174913767L;
	
	private final int INPUT = 0;
	private final int OUTPUT = 1;
//	private HorizontalLayout hlHeader;
	private ToolEditorUI root;
	
	public ToolEditor(ToolEditorUI root) {
		this.root = root;
//		buttonHeader();
		init();
	}
	
//	private void buttonHeader() {
//		hlHeader = new HorizontalLayout();
//		hlHeader.setImmediate(true);
//		Button btAddTool = new Button("Add Tool");
//		btAddTool.setVisible(false);
//		Button btAddInput = new Button("Add Input");
//		Button btAddOutput = new Button("Add Output");
//		Button btAddParameter = new Button("Add Parameter");
//		hlHeader.setSpacing(true);
//		hlHeader.setMargin(true);
//		hlHeader.addComponent(btAddTool);
//		hlHeader.addComponent(btAddInput);
//		hlHeader.addComponent(btAddOutput);
//		hlHeader.addComponent(btAddParameter);
//		this.addComponent(hlHeader);
//		btAddTool.addClickListener(new ClickListener() {
//			private static final long serialVersionUID = 7392773905732834881L;
//
//			@Override
//			public void buttonClick(ClickEvent event) {
//				ToolEditor.this.addTool();
//			}
//		});
//		
//		btAddInput.addClickListener(new ClickListener() {
//			private static final long serialVersionUID = 4949254611179248011L;
//
//			@Override
//			public void buttonClick(ClickEvent event) {
//				ToolEditor.this.addEmptyInput();
//			}
//		});
//		
//		btAddOutput.addClickListener(new ClickListener() {
//			private static final long serialVersionUID = 8773521771491048374L;
//
//			@Override
//			public void buttonClick(ClickEvent event) {
//				ToolEditor.this.addEmptyOutput();
//			}
//		});
//		
//		btAddParameter.addClickListener(new ClickListener() {
//			private static final long serialVersionUID = 4668120288792175832L;
//
//			@Override
//			public void buttonClick(ClickEvent event) {
//				ToolEditor.this.addEmptyParameter();
//			}
//		});
//		
//	}
	
	private void init() {
		this.setImmediate(true);
		addTool();
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
			this.addComponent(new Tool(this));
	}
	
	public void addTool(SADLDescription sadlDescription) {
		Tool tool = new Tool(this).createUIwithData(sadlDescription);
		root.getTreeToolEditor().addTool(tool);
		this.addComponent(tool);
	}
	
	public void addEmptyInput() {
		addEmptyInput(null);
	}
	
	public void addEmptyInput(Input input) {
		this.addComponent(input == null ? new Input(this) : input, getLastIndex(INPUT)+1);
	}
	
	public void addInput(SADLDescription.Input input) {
		Input element = new Input(this);
		root.getTreeToolEditor().addInput(element);
		this.addComponent(element);
		element.fillWithData(input);
		root.getTreeToolEditor().setItemCaption(element, element.toString());
	}
	
	public void addEmptyOutput() {
		addEmptyOutput(null);
	}
	
	public void addEmptyOutput(Output output) {
		this.addComponent(output == null ? new Output(this) : output, getLastIndex(OUTPUT)+1);
	}
	
	public void addOutput(SADLDescription.Output output) {
		Output element = new Output(this);
		root.getTreeToolEditor().addOutput(element);
		this.addComponent(element);
		element.fillWithData(output);
	}
	
	public void addEmptyParameter() {
		addEmptyParameter(null);
	}
	
	public void addEmptyParameter(Parameter parameter) {
		this.addComponent(parameter == null ? new Parameter(this) : parameter);
	}
	
	public void addParameter(SADLDescription.Parameter parameter) {
		Parameter element = new Parameter(this);
		root.getTreeToolEditor().addParameter(element);
		this.addComponent(element);
		element.fillWithData(parameter);
	}
	
	public void removeItems() {
		this.removeAllComponents();
//		this.addComponent(hlHeader);
	}

	public Tool getTool() {
		return root.getTreeToolEditor().getTool();
//		for(int i = 0; i < this.getComponentCount(); i++) {
//			if(this.getComponent(i) instanceof Tool) {
//				return (Tool) this.getComponent(i);
//			}
//		}
//		return null;
	}
	
	public List<SADLDescription.Input> getDaslInputs() {
		ArrayList<SADLDescription.Input> inputs = new ArrayList<SADLDescription.Input>();
		Collection<Input> in = root.getTreeToolEditor().getInputs();
		if(in == null)
			return inputs;
		for(Input input : in) {
			inputs.add(input.getSadlInput());
		}
//		for(int i = 0; i < this.getComponentCount(); i++) {
//			if(this.getComponent(i) instanceof Input) {
//				inputs.add(((Input) this.getComponent(i)).getSadlInput());
//			} else {
//				if(!inputs.isEmpty()) {
//					break;
//				}
//			}
//		}
		return inputs;
	}
	
	public List<SADLDescription.Output> getDaslOutputs() {
		ArrayList<SADLDescription.Output> outputs = new ArrayList<SADLDescription.Output>();
		Collection<Output> out = root.getTreeToolEditor().getOutputs();
		if(out == null)
			return outputs;
		for(Output output : out) {
			outputs.add(output.getSadlOutput());
		}
//		for(int i = 0; i < this.getComponentCount(); i++) {
//			if(this.getComponent(i) instanceof Output) {
//				outputs.add(((Output) this.getComponent(i)).getSadlOutput());
//			} else {
//				if(!outputs.isEmpty()) {
//					break;
//				}
//			}
//		}
		return outputs;
	}
	
	public void addParameters(SADLDescription sadlDescription) {
//		boolean wasParameter = false;
		if(root.getTreeToolEditor().getParameters() == null)
			return;
		for(Parameter parameter : root.getTreeToolEditor().getParameters()) {
			sadlDescription.addParameter(parameter.getSadlParameter());
		}
//		for(int i = 0; i < this.getComponentCount(); i++) {
//			if(this.getComponent(i) instanceof Parameter) {
//				sadlDescription.addParameter(((Parameter) this.getComponent(i)).getSadlParameter());
//				wasParameter = true;
//			} else {
//				if(wasParameter) {
//					break;
//				}
//			}
//		}
//		return sadlDescription;
	}
	
	public void setHeaderToTextEditor(String text) {
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
		if(index+1 < this.getComponentCount()) {
			Component com = this.getComponent(index+1);
			if((component instanceof Input && com instanceof Input) || 
					(component instanceof Output && com instanceof Output) || 
					(component instanceof Parameter && com instanceof Parameter)) {
				this.removeComponent(component);
				this.addComponent(component, index+1);
			}
		}
	}
	
	public ToolEditorUI getToolEditorUI() {
		return root;
	}
}
