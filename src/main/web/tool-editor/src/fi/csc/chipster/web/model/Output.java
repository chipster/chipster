package fi.csc.chipster.web.model;

import com.vaadin.ui.Label;

import fi.csc.chipster.web.tooledit.ToolEditor;
import fi.csc.microarray.description.SADLDescription;

/**
 * Tool editor output model.
 * @author Gintare Pacauskaite
 *
 */
public class Output extends InputOutputUI{

	private static final long serialVersionUID = 6280877684047822425L;
	
	public Output(ToolEditor root) {
		this.root = root;
		createUI();
	}
	
	@Override
	public Output createUI() {
		generateHeader();
		initElements();
		createBasicUI();
		// just for now
		removeTypeRow();
		return this;
	}
	
	@Override
	protected void initElements() {
		super.initElements();
		lbId.setValue("Output file:");
		lbType = new Label("Output is:");
		lbType2 = new Label("Type:");
		
	}
	
	public Output createUIWithData(SADLDescription.Output output) {
		fillWithData(output);
		return this;
	}
	
	/**
	 * Fills up output fields from SADL output
	 * @param output SADL output
	 */
	public void fillWithData(SADLDescription.Output output) {
		name.setValue(output.getName().getDisplayName());
		cbMeta.setValue(output.isMeta());
		description.setValue(getValue(output.getComment()));
		optional.setValue(output.isOptional());
		
		if(output.getName().getPrefix() == null || output.getName().getPrefix().isEmpty()) {
			type.select(SINGLE_FILE);
			id.setValue(output.getName().getID());
			getSingleFileUI();
		} else {
			type.select(MULTI_FILE);
			prefix.setValue(output.getName().getPrefix());
			postfix.setValue(output.getName().getPostfix());
			getMultipleFilesUI();
		}
	}
	
	/**
	 * Creates SADL output from tool editor output
	 * @return SADL output
	 */
	public SADLDescription.Output getSadlOutput() {
		SADLDescription.Output output = new SADLDescription.Output();
		output.setName(getNameFromUI(type.getValue().toString()));
		output.setOptional(optional.getValue());
		output.setComment(getValueOrNull(description.getValue()));
		output.setMeta(cbMeta.getValue());
		
		return output;
	}

	@Override
	protected void generateHeader() {
		lbTitle.setValue(getBoldText("Output"));	
	}
	
	@Override
	public String toString() {
		return getValue((id.getValue() != null && !id.getValue().isEmpty() ? id.getValue() :name.getValue()));
	}
	
	private void removeTypeRow() {
		this.removeRow(this.getComponentArea(type2).getRow1());
	}

}
