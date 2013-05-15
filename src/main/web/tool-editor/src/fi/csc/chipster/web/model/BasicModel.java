package fi.csc.chipster.web.model;

import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextArea;
import com.vaadin.ui.TextField;

public class BasicModel {

	protected TextField name;
	protected TextField id;
	protected TextArea description;
	protected GridLayout grid;
	
	public BasicModel() {
		grid = new GridLayout(3, 4);
		grid.setImmediate(true);
		grid.setMargin(true);
		grid.setSpacing(true);
		grid.setColumnExpandRatio(0, 2);
		grid.setColumnExpandRatio(1, 2);
//		grid.setWidth("100%");
//		grid.setHeight("100%");
		initElements();
		grid.addComponent(name, 0, 0);
		
		grid.addComponent(id, 1,0);
		grid.addComponent(description, 2, 0, 2, 1);
	}
	
	private void initElements() {
		name = new TextField("Display name");
		name.setWidth("100px");
		
		id = new TextField("Id");
		id.setWidth("100px");
		
		description = new TextArea();
		description.setCaption("Description");
		description.setWidth("200px");
	}

}
