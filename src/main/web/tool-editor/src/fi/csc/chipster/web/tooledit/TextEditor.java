package fi.csc.chipster.web.tooledit;

import java.util.List;

import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.TextArea;
import com.vaadin.ui.UI;
import com.vaadin.ui.VerticalLayout;

import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Input;
import fi.csc.microarray.description.SADLDescription.Output;
import fi.csc.microarray.description.SADLDescription.Parameter;
import fi.csc.microarray.description.SADLParser;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class TextEditor extends VerticalLayout{
	private static final long serialVersionUID = -7074541336842177583L;
	
	private TextArea txtArea;
	private Button btUpdate;
	private ToolEditorUI root;
	
	public static final String NEW_LINE = "\n";
	

	public TextEditor(ToolEditorUI root) {
		this.root = root;
		init();
	}
	
	private void init() {
		
		btUpdate = new Button("Update");
		btUpdate.addClickListener(new ClickListener() {
			private static final long serialVersionUID = 1487893808578560989L;

			@Override
			public void buttonClick(ClickEvent event) {
				String text = takeHeader(getText());
				
				ChipsterSADLParser parser = new ChipsterSADLParser();
				try {
					SADLDescription description = parser.parse(text);
//					System.out.println(description);
					root.getToolEditor().removeItems();
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
					
				} catch (ParseException e) {
					e.printStackTrace();
				}
				
			}
		});
		this.addComponent(btUpdate);
		
		txtArea = new TextArea();
		txtArea.setRows(20);
		txtArea.setWidth("80%");
		
		this.addComponent(txtArea);
	}
	public void setText(String text) {
		txtArea.setValue(replaceHeader(text));
	}
	
	public String getText() {
		return txtArea.getValue();
	}
	
	public String takeHeader(String input) {
		StringBuilder output = new StringBuilder();
		input = input.replace("#", "");
		String[] array = input.split(NEW_LINE);
		int i = 0;
		array[i] = array[i].trim();
//		System.out.println(array.length + " " + array[0]);
		while (i < array.length && (array[i].startsWith("TOOL") || array[i].startsWith("INPUT") 
				|| array[i].startsWith("OUTPUT") || array[i].startsWith("PARAMETER"))) {
			output.append(array[i] + NEW_LINE);
			i++;
			if(i < array.length)
				array[i] = array[i].trim();
		}
//		System.out.println("output:" + output);
		return output.toString();
	}
	
	private String replaceHeader(String header) {
		StringBuilder old = new StringBuilder();
		String text = txtArea.getValue();
		String[] array = text.split(NEW_LINE);
		int i = 0;
		while (i < array.length && (array[i].startsWith("# TOOL") || array[i].startsWith("# INPUT") 
				|| array[i].startsWith("# OUTPUT") || array[i].startsWith("# PARAMETER"))) {
			old.append(array[i] + NEW_LINE);
			i++;
		}
		return text.replace(old, header);
	}

}
