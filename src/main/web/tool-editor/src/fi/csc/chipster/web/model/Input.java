package fi.csc.chipster.web.model;

import com.vaadin.ui.Label;

import fi.csc.chipster.web.listener.CSCTypeValueChangeListener;
import fi.csc.chipster.web.tooledit.ToolEditor;
import fi.csc.microarray.description.GenericInputTypes;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;

/**
 * Input for tool editor
 * @author Gintare Pacauskaite
 *
 */
public class Input extends InputOutputUI{
	
	private static final long serialVersionUID = -7370540274349829538L;
	
	public Input(ToolEditor root) {
		this.root = root;
		createUI();
	}
	
	@Override
	public Input createUI() {
		generateHeader();
		initElements();
		createBasicUI();
		return this;
	}
	
	@Override
	protected void generateHeader() {
		lbTitle.setValue(getBoldText("Input:"));
	}
	
	@Override
	protected void initElements() {
		super.initElements();
		
		lbId.setValue("Input file:");
		lbType = new Label("Input is:");
		lbType2 = new Label("Type:");
		addInputTypes();
		type2.setNullSelectionAllowed(false);
		type2.select(GenericInputTypes.GENERIC);
		type2.addValueChangeListener(new CSCTypeValueChangeListener(this));
	}
	
	public Input createUIWithData(SADLDescription.Input input) {
		fillWithData(input);
		return this;
	}
	
	/**
	 * Fills input fields with data from SADL input.
	 * @param input
	 */
	public void fillWithData(SADLDescription.Input input) {
		type2.select(input.getType());
		name.setValue(input.getName().getDisplayName());
		description.setValue(getValue(input.getComment()));
		cbMeta.setValue(input.isMeta());
		optional.setValue(input.isOptional());
		
		if(input.getName().getPrefix() == null || input.getName().getPrefix().isEmpty()) {
			type.select(SINGLE_FILE);
			id.setValue(input.getName().getID());
			getSingleFileUI();
		} else {
			type.select(MULTI_FILE);
			prefix.setValue(input.getName().getPrefix());
			postfix.setValue(input.getName().getPostfix());
			getMultipleFilesUI();
		}
	}
	
	/**
	 * Creates SADL input from tool editor input fields
	 * @return SADL input
	 */
	public SADLDescription.Input getSadlInput() {
		SADLDescription.Input input = new SADLDescription.Input();
		input.setName(getNameFromUI(type.getValue().toString()));
		input.setOptional(optional.getValue());
		input.setComment(getValueOrNull(description.getValue()));
		input.setType((InputType) type2.getValue());
		input.setMeta(cbMeta.getValue());
		return input;
	}
	
	private void addInputTypes() {
		type2.addItem(ChipsterInputTypes.AFFY);
		type2.setItemCaption(ChipsterInputTypes.AFFY, ChipsterInputTypes.AFFY.getName());
		type2.addItem(ChipsterInputTypes.BAM);
		type2.setItemCaption(ChipsterInputTypes.BAM, ChipsterInputTypes.BAM.getName());
		type2.addItem(ChipsterInputTypes.CDNA);
		type2.setItemCaption(ChipsterInputTypes.CDNA, ChipsterInputTypes.CDNA.getName());
		type2.addItem(ChipsterInputTypes.FASTA);
		type2.setItemCaption(ChipsterInputTypes.FASTA, ChipsterInputTypes.FASTA.getName());
		type2.addItem(ChipsterInputTypes.GENE_EXPRS);
		type2.setItemCaption(ChipsterInputTypes.GENE_EXPRS, ChipsterInputTypes.GENE_EXPRS.getName());
		type2.addItem(ChipsterInputTypes.GENELIST);
		type2.setItemCaption(ChipsterInputTypes.GENELIST, ChipsterInputTypes.GENELIST.getName());
		type2.addItem(ChipsterInputTypes.GTF);
		type2.setItemCaption(ChipsterInputTypes.GTF, ChipsterInputTypes.GTF.getName());
		type2.addItem(ChipsterInputTypes.PHENODATA);
		type2.setItemCaption(ChipsterInputTypes.PHENODATA, ChipsterInputTypes.PHENODATA.getName());
		type2.addItem(GenericInputTypes.GENERIC);
		type2.setItemCaption(GenericInputTypes.GENERIC, GenericInputTypes.GENERIC.getName());
	}
	
	@Override
	public String toString() {
		return getValue((id.getValue() != null && !id.getValue().isEmpty() ? id.getValue() :name.getValue())) + " " + getValue(getType());
	}

	private boolean isMultipleInput() {
		if (MULTI_FILE.equals(type.getValue().toString())) {
			return true;
		}
		return false;
	}
	
	@Override
	public boolean isValid() {
		if (isMultipleInput()) {
			if (!prefix.getValue().isEmpty() && !postfix.getValue().isEmpty()) {
				return true;
			} else {
				return false;
			}
		} else {
			return super.isValid();
		}
	}

}
