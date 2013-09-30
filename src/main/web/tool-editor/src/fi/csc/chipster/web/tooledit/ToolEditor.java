package fi.csc.chipster.web.tooledit;

import java.io.IOException;
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

/**
 * Tool editor class.
 * 
 * @author Gintare Pacauskaite
 * 
 */
public class ToolEditor extends VerticalLayout {

	private static final long serialVersionUID = -4777411581174913767L;

	private final int INPUT = 0;
	private final int OUTPUT = 1;
	private ToolEditorUI root;

	public ToolEditor(ToolEditorUI root) {
		this.root = root;
		init();
		root.getTreeToolEditor().addTool(new Tool(this));
	}

	private void init() {
		this.setImmediate(true);
		addTool();
	}

	public void addTool() {
		if (getTool() == null)
			this.addComponent(new Tool(this));
	}

	/**
	 * Adds tool from SADL Description to tool editor and summary tree
	 * 
	 * @param sadlDescription
	 */
	public void addTool(SADLDescription sadlDescription) {
		Tool tool = new Tool(this).createUIwithData(sadlDescription);
		tool.setTitleDescriptionValue(tool.getTitle());
		root.getTreeToolEditor().addTool(tool);
		root.getTreeToolEditor().updateToolTitle(tool.getTitle());
		this.addComponent(tool);
	}

	public void addEmptyInput() {
		addEmptyInput(null);
	}

	/**
	 * Adds new input in the rigth place
	 * 
	 * @param input
	 */
	public void addEmptyInput(Input input) {
		this.addComponent(input == null ? new Input(this) : input,
				getLastIndex(INPUT) + 1);
	}

	/**
	 * Adds SADL input to tool editor and summary tree
	 * 
	 * @param input
	 */
	public void addInput(SADLDescription.Input input) {
		Input element = new Input(this);
		root.getTreeToolEditor().addInput(element);
		this.addComponent(element);
		element.fillWithData(input);
		element.setTitleDescriptionValue(element.toString());
		root.getTreeToolEditor().setItemCaption(element, element.toString());
	}

	/**
	 * Adds empty output
	 */
	public void addEmptyOutput() {
		addEmptyOutput(null);
	}

	/**
	 * Adds empty output in the rigth place
	 * 
	 * @param output
	 */
	public void addEmptyOutput(Output output) {
		this.addComponent(output == null ? new Output(this) : output,
				getLastIndex(OUTPUT) + 1);
	}

	/**
	 * Adds SADL output to tool editor and summary tree
	 * 
	 * @param output
	 */
	public void addOutput(SADLDescription.Output output) {
		Output element = new Output(this);
		root.getTreeToolEditor().addOutput(element);
		this.addComponent(element);
		element.fillWithData(output);
		element.setTitleDescriptionValue(element.toString());
		root.getTreeToolEditor().setItemCaption(element, element.toString());
	}

	/**
	 * Adds empty parameter to tool editor
	 */
	public void addEmptyParameter() {
		addEmptyParameter(null);
	}

	/**
	 * Adds parameter to the rigth place. If parameter is null, adds new
	 * parameter to tool editor
	 * 
	 * @param parameter
	 */
	public void addEmptyParameter(Parameter parameter) {
		this.addComponent(parameter == null ? new Parameter(this) : parameter);
	}

	/**
	 * Adds SADL parameter to tool editor and summary
	 * 
	 * @param parameter
	 */
	public void addParameter(SADLDescription.Parameter parameter) {
		Parameter element = new Parameter(this);
		root.getTreeToolEditor().addParameter(element);
		this.addComponent(element);
		element.fillWithData(parameter);
		element.setTitleDescriptionValue(element.toString());
		root.getTreeToolEditor().setItemCaption(element, element.toString());
	}

	/**
	 * Removes all element from tool editor
	 */
	public void removeItems() {
		this.removeAllComponents();
	}

	public Tool getTool() {
		return root.getTreeToolEditor().getTool();
	}

	/**
	 * Generates SADL inputs from tool editor inputs
	 * 
	 * @return list of SADL inputs, or empty list if no inputs were found in
	 *         summary tree
	 * @throws Exception
	 *             if input is not valid
	 */
	public List<SADLDescription.Input> getSADLInputs() throws Exception {
		ArrayList<SADLDescription.Input> inputs = new ArrayList<SADLDescription.Input>();
		Collection<Input> in = root.getTreeToolEditor().getInputs();
		if (in == null)
			return inputs;
		for (Input input : in) {
			if (input.isValid())
				inputs.add(input.getSadlInput());
			else
				throw new Exception("Invalid input: " + input);
		}
		return inputs;
	}

	/**
	 * Generates SADL outputs from tool editor outputs.
	 * 
	 * @return list of SADL outputs, or empty list if no outputs were found in
	 *         summary tree.
	 * @throws Exception
	 *             output is valid
	 */
	public List<SADLDescription.Output> getSADLOutputs() throws Exception {
		ArrayList<SADLDescription.Output> outputs = new ArrayList<SADLDescription.Output>();
		Collection<Output> out = root.getTreeToolEditor().getOutputs();
		if (out == null)
			return outputs;
		for (Output output : out) {
			if (output.isValid())
				outputs.add(output.getSadlOutput());
			else
				throw new Exception();
		}
		return outputs;
	}

	/**
	 * Generates and add SADL parameters to SADL description from tool editor
	 * parameters
	 * 
	 * @param sadlDescription
	 *            SADL description with parameters
	 * @throws Exception
	 *             if parameter is not valid
	 */
	public void addParameters(SADLDescription sadlDescription) throws Exception {
		if (root.getTreeToolEditor().getParameters() == null)
			return;
		for (Parameter parameter : root.getTreeToolEditor().getParameters()) {
			if (parameter.isValid())
				sadlDescription.addParameter(parameter.getSadlParameter());
			else
				throw new Exception();
		}
	}

	/**
	 * Sets header for text editor
	 * 
	 * @param text header without prefix
	 * @throws IOException 
	 */
	public void setHeaderToTextEditor(String text) throws IOException {
		root.getTextEditor().setHeader(text);
	}


	private int getLastIndex(int type) {
		int toolIndex = -1;
		int outputIndex = -1;
		int inputIndex = -1;
		for (int i = 0; i < this.getComponentCount(); i++) {
			if (this.getComponent(i) instanceof Input) {
				inputIndex = i;
			} else if (this.getComponent(i) instanceof Output) {
				outputIndex = i;
			} else if (this.getComponent(i) instanceof Tool) {
				toolIndex = i;
			} else if (this.getComponent(i) instanceof Parameter) {
				break;
			}
		}
		if (type == INPUT) {
			if (inputIndex > -1)
				return inputIndex;
			else if (toolIndex > -1)
				return toolIndex;
			else
				return 0;
		} else if (type == OUTPUT) {
			if (outputIndex > -1)
				return outputIndex;
			else if (inputIndex > -1)
				return inputIndex;
			else if (toolIndex > -1)
				return toolIndex;
			else
				return 0;
		} else
			return 0;
	}

	@Deprecated
	public void moveUpComponent(BasicModel component) {
		int index = this.getComponentIndex(component);
		if (index - 1 > 0) {
			Component com = this.getComponent(index - 1);
			if ((component instanceof Input && com instanceof Input)
					|| (component instanceof Output && com instanceof Output)
					|| (component instanceof Parameter && com instanceof Parameter)) {
				this.removeComponent(component);
				this.addComponent(component, index - 1);
			}
		}
	}

	@Deprecated
	public void moveDownComponent(BasicModel component) {
		int index = this.getComponentIndex(component);
		if (index + 1 < this.getComponentCount()) {
			Component com = this.getComponent(index + 1);
			if ((component instanceof Input && com instanceof Input)
					|| (component instanceof Output && com instanceof Output)
					|| (component instanceof Parameter && com instanceof Parameter)) {
				this.removeComponent(component);
				this.addComponent(component, index + 1);
			}
		}
	}

	public ToolEditorUI getToolEditorUI() {
		return root;
	}
}
