package fi.csc.chipster.web.tooledit;

import com.vaadin.ui.TextArea;
import com.vaadin.ui.VerticalLayout;

public class TextEditor extends VerticalLayout{
	private static final long serialVersionUID = -7074541336842177583L;
	
	private TextArea txtArea;
	

	public TextEditor() {
		init();
	}
	
	private void init() {
		txtArea = new TextArea();
		txtArea.setRows(20);
		txtArea.setWidth("80%");
		
		this.addComponent(txtArea);
	}

}
