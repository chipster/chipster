package fi.csc.chipster.web.model;

import com.vaadin.ui.ComboBox;
import com.vaadin.ui.Label;

import fi.csc.chipster.web.tooledit.ToolEditor;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLDescription.Name;

public class Tool extends BasicModel{
	
	private static final long serialVersionUID = -6173788013713655305L;
	private Label lbModule;
	private Label lbCategory;
	
	private ComboBox module;
	private ComboBox category;
	
	public Tool(ToolEditor root) {
		this.root = root;
		createUI();
	}
	
	@Override
	public Tool createUI() {
		generateHeader();
		initElements();
		addRow(lbModule, module);
		addRow(lbCategory, category);
		addRow(lbDescription, description);
		return this;
	}
	
	@Override
	protected void generateHeader() {
		lbTitle.setValue(getBoldText("Tool"));
		initFooter();
		btUp.setVisible(false);
		btDown.setVisible(false);
		btDelete.setVisible(false);
		this.replaceComponent(lbTitleDescription, hLayoutTitle);
	}
	
	public Tool createUIwithData(SADLDescription sadlDescription) {
//		createUI();
		fillWithData(sadlDescription);
		return this;
	}
	
	private void initElements() {
		lbId.setValue("Tool name:");
		lbModule = new Label("Module:");
		lbCategory = new Label("Category:");
		
		module = new ComboBox();
		module.setWidth(WIDTH);
		
		category = new ComboBox();
		category.setWidth(WIDTH);
		
		id.setRequired(true);
		name.setRequired(true);
		description.setRequired(true);
	}
	
	
	private void fillWithData(SADLDescription sadlDescription) {
		id.setValue(sadlDescription.getName().getID());
		name.setValue(sadlDescription.getName().getDisplayName());
		description.setValue(getValue(sadlDescription.getComment()));
	}
	
	public SADLDescription createSadlDescription() {
		if(id.isValid() && name.isValid() && description.isValid()) {
			SADLDescription sadlDescription = new SADLDescription(Name.createName(id.getValue(), name.getValue()));
			sadlDescription.setComment(getValueOrNull(description.getValue()));
			return sadlDescription;
		} else {
			return null;
		}
	}

	@Override
	protected String getType() {
		// TODO Auto-generated method stub
		return null;
	}
}
