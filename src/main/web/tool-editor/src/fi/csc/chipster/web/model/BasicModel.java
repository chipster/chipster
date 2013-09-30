package fi.csc.chipster.web.model;

import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Alignment;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.Component;
import com.vaadin.ui.GridLayout;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TextArea;
import com.vaadin.ui.TextField;

import fi.csc.chipster.web.listener.CSCTextChangeListener;
import fi.csc.chipster.web.tooledit.ToolEditor;
import fi.csc.microarray.description.SADLDescription.Name;

/**
 * Basic model for tool, input, output, parameter
 * @author Gintare Pacauskaite
 *
 */
public abstract class BasicModel extends GridLayout{
	
	private static final long serialVersionUID = 7942793110778486624L;
	public static final String WIDTH = "300px";
	public static final String OPTIONAL = "optional";
	public static final String NOT_OPTIONAL = "not optional";
	public static final String SINGLE_FILE = "Single file";
	public static final String MULTI_FILE = "Multiple files";
	public static final String MULTI_FILE_TEXT = "{...}";
	public static final String REQUIRED_TEXT = "this field cannot be empty";
	
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
		name.setDescription("Display name for the element");
		name.setImmediate(true);
		name.addTextChangeListener(new CSCTextChangeListener(this));
		
		id = new TextField();
		id.setDescription("file name or unique identification");
		id.setImmediate(true);
		id.setRequired(true);
		id.setRequiredError(REQUIRED_TEXT);
		id.setWidth(WIDTH);
		id.addTextChangeListener(new CSCTextChangeListener(this, true));
		
		description = new TextArea();
		description.setWidth(WIDTH);
		description.setDescription("Short description");
		
		layout = new HorizontalLayout();
		
		prefix = new TextField();
		postfix = new TextField();
		
		layout.addComponent(prefix);
		layout.addComponent(new Label(MULTI_FILE_TEXT));
		layout.addComponent(postfix);
		
		optional = new CheckBox();
		optional.setDescription("Is this element optional");
		optional.setImmediate(true);
	}
	
	abstract protected String getType();
	
	public String getTypeValue() {
		return getType();
	}
	
	public String getNameValue() {
		return name.getValue().toString();
	}
	
	/**
	 * sets title
	 * @param text 
	 */
	public void setTitleDescriptionValue(String text) {
		lbTitleDescription.setValue(getBoldText(text));
	}
	
	/**
	 * adds row
	 * @param label
	 * @param component
	 */
	protected void addRow(Component label, Component component) {
		this.addComponent(label);
		this.addComponent(component);
		this.setComponentAlignment(label, Alignment.MIDDLE_LEFT);
	}
	
	/**
	 * get string value
	 * @param value
	 * @return value or empty string
	 */
	protected String getValue(String value) {
		return value == null || value.isEmpty() ? "" : value;
	}
	
	/**
	 *  get string value
	 * @param value
	 * @return value or null
	 */
	protected String getValueOrNull(String value) {
		return value == null || value.isEmpty() ? null : value;
	}
	
	public abstract GridLayout createUI();
	protected abstract void generateHeader();
	
	/**
	 * Creates UI for multiple files
	 */
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
	
	/**
	 * Gets name from UI
	 * @param value
	 * @return name, or null if some other type than sigle file / multipe file was selected
	 */
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
	
	public String getId() {
		return id.getValue();
	}
	
	/**
	 * Input overrides this to avoid multiple input being invalid
	 * 
	 * @return
	 */
	public boolean isValid() {
		return id.isValid();
	}

}
