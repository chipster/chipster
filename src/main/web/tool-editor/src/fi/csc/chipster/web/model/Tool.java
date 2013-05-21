package fi.csc.chipster.web.model;

import com.vaadin.ui.ComboBox;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.Label;

import fi.csc.microarray.description.SADLDescription;

public class Tool extends BasicModel{

//	private String name;
//	private String id;
//	private String description;
	
	private Label lbModule;
	private Label lbCategory;
	
	private ComboBox module;
	private ComboBox category;
	
	@Override
	public GridLayout createUI() {
		
//		GridLayout grid = new GridLayout(3, 3);
//		grid.setImmediate(true);
//		grid.setMargin(true);
//		grid.setSpacing(true);
//		grid.setColumnExpandRatio(0, 2);
//		grid.setColumnExpandRatio(1, 2);
//		grid.setWidth("100%");
//		grid.setHeight("100%");
//		grid.addComponent(new Label("Tool"), 0, 0);
		initElements();
		addRow(lbModule, module);
		addRow(lbCategory, category);
		addRow(lbDescription, description);
		
//		grid.addComponent(name, 0, 1);
//		
//		grid.addComponent(id, 1,1);
//		grid.addComponent(description, 2, 1, 2, 2);
//		grid.addComponent(module, 0, 1);
//		grid.addComponent(category, 1, 1);
//		VerticalLayout vertical = new VerticalLayout();
////		vertical.setWidth("100%");
//		vertical.addComponent(new Label("Tool"));
//		HorizontalLayout horizontal = new HorizontalLayout();
//		horizontal.setMargin(true);
//		horizontal.setSpacing(true);
//		initElements();
//		horizontal.addComponent(name);
//		horizontal.addComponent(id);
//		horizontal.addComponent(description);
//		vertical.addComponent(horizontal);
		return grid;
	}
	
	public GridLayout createUIwithData(SADLDescription sadlDescription) {
		createUI();
		fillWithData(sadlDescription);
		return grid;
	}
	
	private void initElements() {
//		name = new TextField("Display name");
//		name.setWidth("100%");
//		
//		id = new TextField("Id");
//		id.setWidth("100%");
//		
//		description = new TextArea();
//		description.setCaption("Description");
//		description.setWidth("200px");
		lbId.setValue("Tool name:");
		lbModule = new Label("Module:");
		lbCategory = new Label("Category:");
		
		module = new ComboBox();
		module.setWidth("100%");
		
		category = new ComboBox();
		category.setWidth("100%");
	}
	
	
	private void fillWithData(SADLDescription sadlDescription) {
		id.setValue(sadlDescription.getName().getID());
		name.setValue(sadlDescription.getName().getDisplayName());
		description.setValue(sadlDescription.getComment());
	}
}
