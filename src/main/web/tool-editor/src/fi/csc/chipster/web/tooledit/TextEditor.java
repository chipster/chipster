package fi.csc.chipster.web.tooledit;

import com.vaadin.ui.Alignment;
import com.vaadin.ui.Button;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.TextArea;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.listener.CSCTextToToolClickListener;
import fi.csc.chipster.web.listener.CSCToolToTextClickListener;

public class TextEditor extends VerticalLayout{
	private static final long serialVersionUID = -7074541336842177583L;
	
	private TextArea txtArea;
	private Button btUpdateToolEditor;
	private Button btUpdateTextEditor;
	private ToolEditorUI root;
	
	public static final String NEW_LINE = "\n";
	

	public TextEditor(ToolEditorUI root) {
		this.root = root;
		init();
	}
	
	private void init() {
		
		HorizontalLayout hLayout = new HorizontalLayout();
		hLayout.setSpacing(true);
		
		btUpdateTextEditor = new Button();
		btUpdateTextEditor.setIcon(Icon.getResource(Icon.getDownButtonIconPath()));
		btUpdateTextEditor.addClickListener(new CSCToolToTextClickListener(root));
		hLayout.addComponent(btUpdateTextEditor);
		btUpdateToolEditor = new Button();
		btUpdateToolEditor.setIcon(Icon.getResource(Icon.getUpButtonIconPath()));
		btUpdateToolEditor.addClickListener(new CSCTextToToolClickListener(root));
		hLayout.addComponent(btUpdateToolEditor);
		this.addComponent(hLayout);
		this.setComponentAlignment(hLayout, Alignment.MIDDLE_CENTER);
		
		txtArea = new TextArea();
		txtArea.setRows(20);
		txtArea.setWidth("100%");
		
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
		StringBuilder newText = new StringBuilder();
		String text = txtArea.getValue();
		String[] array = text.split(NEW_LINE);
		int i = 0;
		while (i < array.length && (array[i].startsWith("# TOOL") || array[i].startsWith("# INPUT") 
				|| array[i].startsWith("# OUTPUT") || array[i].startsWith("# PARAMETER"))) {
			i++;
		}
		newText.append(header);
		for(int index = i; index < array.length; index++) {
			newText.append(array[index] + NEW_LINE);
		}
		return newText.toString();
	}

}
