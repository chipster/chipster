package fi.csc.chipster.web.model;

import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Alignment;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.Component;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextArea;
import com.vaadin.ui.TextField;

import fi.csc.chipster.web.listener.CSCTextChangeListener;
import fi.csc.chipster.web.tooledit.Icon;
import fi.csc.chipster.web.tooledit.ToolEditor;
import fi.csc.microarray.description.SADLDescription.Name;

public abstract class BasicModel extends GridLayout{
	
	private static final long serialVersionUID = 7942793110778486624L;
	public static final String WIDTH = "300px";
	public static final String OPTIONAL = "optional";
	public static final String NOT_OPTIONAL = "not optional";
	public static final String SINGLE_FILE = "Single file";
	public static final String MULTI_FILE = "Multiple files";
	public static final String MULTI_FILE_TEXT = "{...}";
	
	protected Label lbName;
	protected Label lbId;
	protected Label lbDescription;
	protected Label lbTitle;
	protected Label lbTitleDescription;
	protected Label lbOptional;
	
	
	protected TextField name;
	protected TextField id;
	protected TextArea description;
	protected TextField prefix;
	protected TextField postfix;
	protected Button btUp;
	protected Button btDown;
	protected Button btDelete;
	protected HorizontalLayout hLayoutTitle;
	protected CheckBox optional;
	
	protected HorizontalLayout layout;
	protected ToolEditor root;
	
	public BasicModel() {
		this.setColumns(2);
		this.setImmediate(true);
		this.setMargin(true);
		this.setSpacing(true);
		this.setColumnExpandRatio(0, 20);
		initHeadeer();
		initElements();
		addRow(lbId, id);
		addRow(lbName, name);
	}
	
	
	protected void initFooter() {
		btUp = new Button();
		btDown = new Button();
		btDelete = new Button();
		btDelete.addClickListener(new ClickListener() {
			private static final long serialVersionUID = -689182485982297845L;

			@Override
			public void buttonClick(ClickEvent event) {
				root.removeComponent(BasicModel.this);
				root.getToolEditorUI().getTreeToolEditor().removeItem(BasicModel.this);
			}
		});
		btUp.addClickListener(new ClickListener() {
			private static final long serialVersionUID = -5641105378128102121L;

			@Override
			public void buttonClick(ClickEvent event) {
				root.moveUpComponent(BasicModel.this);
			}
		});
		btDown.addClickListener(new ClickListener() {
			private static final long serialVersionUID = 5198164039029102629L;

			@Override
			public void buttonClick(ClickEvent event) {
				root.moveDownComponent(BasicModel.this);
			}
		});
		String width = "50px";
		btDelete.setIcon(Icon.getResource(Icon.getDeleteButtonIconPath()));
		btDelete.setImmediate(true);
		btDelete.setWidth(width);
		btUp.setImmediate(true);
		btUp.setIcon(Icon.getResource(Icon.getUpButtonIconPath()));
		btUp.setWidth(width);
		btDown.setImmediate(true);
		btDown.setIcon(Icon.getResource(Icon.getDownButtonIconPath()));
		btDown.setWidth(width);
		hLayoutTitle = new HorizontalLayout();
		hLayoutTitle.setSpacing(true);
//		hLayoutTitle.setMargin(true);
//		hLayoutTitle.addComponent(btUp);
//		hLayoutTitle.addComponent(btDown);
		hLayoutTitle.addComponent(btDelete);
	}
	
	
	protected void generateFooter() {
		initFooter();
		addRow(new Label(), hLayoutTitle);
//		this.addComponent(hLayoutTitle);
	}
	
	
	private void initHeadeer() {
		lbTitle = new Label();
		lbTitle.setContentMode(ContentMode.HTML);
		lbTitleDescription = new Label();
		lbTitleDescription.setContentMode(ContentMode.HTML);
		lbTitleDescription.setImmediate(true);
		addRow(lbTitle, lbTitleDescription);
	}
	
	private void initElements() {
		
		lbName = new Label("Display name:");
		lbId = new Label();
		lbDescription = new Label("Description:");
		lbOptional = new Label("Optional:");
		
		name = new TextField();
		name.setWidth(WIDTH);
		name.setImmediate(true);
		name.addTextChangeListener(new CSCTextChangeListener(this));
		
		id = new TextField();
		id.setWidth(WIDTH);
		
		description = new TextArea();
		description.setWidth(WIDTH);
		
		layout = new HorizontalLayout();
		
		prefix = new TextField();
		postfix = new TextField();
		
		layout.addComponent(prefix);
		layout.addComponent(new Label(MULTI_FILE_TEXT));
		layout.addComponent(postfix);
		
		optional = new CheckBox();
		optional.setImmediate(true);
		optional.addValueChangeListener(new ValueChangeListener() {
			private static final long serialVersionUID = 1L;

			@Override
			public void valueChange(ValueChangeEvent event) {
				root.getToolEditorUI().getTreeToolEditor().setItemCaption(BasicModel.this, BasicModel.this.toString());
			}
		});
	}
	
	abstract protected String getType();
	
	public String getTypeValue() {
		return getType();
	}
	
	public String getNameValue() {
		return name.getValue().toString();
	}
	
	public void setTitleDescriptionValue(String text) {
		lbTitleDescription.setValue(getBoldText(text));
	}
	
	protected void addRow(Component label, Component component) {
		this.addComponent(label);
		this.addComponent(component);
		this.setComponentAlignment(label, Alignment.MIDDLE_LEFT);
	}
	
	protected String getValue(String value) {
		return value == null || value.isEmpty() ? "" : value;
	}
	
	protected String getValueOrNull(String value) {
		return value == null || value.isEmpty() ? null : value;
	}
	
	public abstract GridLayout createUI();
	protected abstract void generateHeader();
	
	protected void getMultipleFilesUI() {
		Area area = this.getComponentArea(id);
		if(area == null)
			return;
		this.removeComponent(id);
		this.addComponent(layout, area.getColumn1(), area.getRow1());
		
	}
	
	protected void getSingleFileUI() {
		Area area = this.getComponentArea(layout);
		if(area == null)
			return;
		this.removeComponent(layout);
		this.addComponent(id, area.getColumn1(), area.getRow1());
	}
	
	protected Name getNameFromUI(String value) {
		if(value.equals(SINGLE_FILE)) {
			return Name.createName(id.getValue(), name.getValue());
		} else if(value.equals(MULTI_FILE)) {
			return Name.createNameSet(prefix.getValue(), postfix.getValue(), name.getValue());
		} else {
			return null;
		}
	}
	
	protected String getBoldText(String text) {
		return "<b>" + text + "</b>";
	}
	
	public ToolEditor getToolEditor() {
		return root;
	}
	
	public boolean isOptional() {
		return optional.getValue();
	}

}
