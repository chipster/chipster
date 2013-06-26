package fi.csc.chipster.web.model;

import com.vaadin.ui.ComboBox;
import com.vaadin.ui.Label;

import fi.csc.chipster.web.listener.CSCTextChangeListener;
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
//		initFooter();
//		btUp.setVisible(false);
//		btDown.setVisible(false);
//		btDelete.setVisible(false);
		//this.replaceComponent(lbTitleDescription, hLayoutTitle);
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
		
		id.setImmediate(true);
		id.addTextChangeListener(new CSCTextChangeListener(this, true));
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
	
	public String getTitle() {
		return getTitle(module.getValue(), category.getValue(), name.getValue(), id.getValue());
	}
	
	public String getTitleId(String id) {
		return getTitle(module.getValue(), category.getValue(), name.getValue(), id);
	}
	
	public String getTitleDisplayName(String displayName) {
		return getTitle(module.getValue(), category.getValue(), displayName, id.getValue());
	}
	
	public String getTitle(Object module, Object category, String displayName, String id) {
		StringBuilder title = new StringBuilder();
		if(module != null) {
			title.append(module);
			title.append(" /");
		}
		if(category != null) {
			title.append(" " + category);
			title.append(" /");
		}
		if(displayName != null && !displayName.isEmpty())
			title.append(" " + displayName);
		if(id != null && !id.isEmpty()) {
			title.append(" (" + id + ")");
		}
		return title.toString();
	}
	
	public void setToTitle(String title) {
		lbTitleDescription.setValue(getBoldText(title));
	}
}
