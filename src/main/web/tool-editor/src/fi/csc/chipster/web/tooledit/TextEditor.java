package fi.csc.chipster.web.tooledit;

import com.vaadin.ui.TextArea;
import com.vaadin.ui.VerticalLayout;

/**
 * Text editor
 * @author Gintare Pacauskaite
 *
 */
public class TextEditor extends VerticalLayout{
	private static final long serialVersionUID = -7074541336842177583L;
	
	private TextArea txtArea;
	
	public static final String NEW_LINE = "\n";

	public TextEditor(ToolEditorUI root) {
		init();
	}
	
	private void init() {
		
		txtArea = new TextArea();
		txtArea.setSizeFull();
		// for some reasons size full does not do anything to height
		txtArea.setRows(50);
		
		this.addComponent(txtArea);
	}
	public void setText(String text) {
		txtArea.setValue(replaceHeader(text));
	}
	
	public String getText() {
		return txtArea.getValue();
	}
	
	/**
	 * Takes header from text
	 * @param input
	 * @return tool header
	 */
	public String takeHeader(String input) {
		if(input == null || input.isEmpty())
			return "";
		StringBuilder output = new StringBuilder();
		input = input.replace("#", "");
		String[] array = input.split(NEW_LINE);
		int i = 0;
		array[i] = array[i].trim();
		while (i < array.length && (array[i].startsWith("TOOL") || array[i].startsWith("INPUT") 
				|| array[i].startsWith("OUTPUT") || array[i].startsWith("PARAMETER"))) {
			output.append(array[i] + NEW_LINE);
			i++;
			if(i < array.length)
				array[i] = array[i].trim();
		}
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
	
	public void clearAllText() {
		txtArea.setValue("");
	}

}
