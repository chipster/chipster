package fi.csc.chipster.web.model;

import com.vaadin.ui.Component;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextArea;
import com.vaadin.ui.TextField;

public abstract class BasicModel {
	
	public static final String OPTIONAL = "optional";
	public static final String NOT_OPTIONAL = "not optional";
	
	protected Label lbName;
	protected Label lbId;
	protected Label lbDescription;
	
	protected TextField name;
	protected TextField id;
	protected TextArea description;
	
	protected GridLayout grid;
	
	public BasicModel() {
		grid = new GridLayout();
		grid.setColumns(2);
		grid.setImmediate(true);
		grid.setMargin(true);
		grid.setSpacing(true);
		grid.setColumnExpandRatio(0, 2);
//		grid.setColumnExpandRatio(1, 2);
//		grid.setWidth("100%");
//		grid.setHeight("100%");
		initElements();
		addRow(lbId, id);
		addRow(lbName, name);
//		grid.addComponent(name, 0, 0);
//		
//		grid.addComponent(id, 1,0);
//		grid.addComponent(description, 2, 0, 2, 1);
	}
	
	private void initElements() {
		
		lbName = new Label("Display name:");
		lbId = new Label();
		lbDescription = new Label("Description");
		
		name = new TextField();
		name.setWidth("200px");
		
		id = new TextField();
		id.setWidth("200px");
		
		description = new TextArea();
		description.setWidth("200px");
	}
	
	protected void addRow(Label label, Component component) {
		grid.addComponent(label);
		grid.addComponent(component);
	}
	
	public abstract GridLayout createUI();

}
