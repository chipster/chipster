package fi.csc.chipster.web.model;

import com.vaadin.ui.Label;

import fi.csc.microarray.description.GenericInputTypes;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLSyntax.InputType;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;

public class Input extends InputOutputUI{
	
	private static final long serialVersionUID = -7370540274349829538L;
	
	
	@Override
	public Input createUI() {
//		grid.addComponent(new Label("Input"), 0, 0);
		initElements();
		createBasicUI();
		return this;
	}
	
	@Override
	protected void initElements() {
		super.initElements();
		
		lbId.setValue("Input file:");
		lbOptional = new Label("Input is:");
		lbType = new Label("Input is:");
		lbType2 = new Label("Type:");
		addInputTypes();
		type2.setNullSelectionAllowed(false);
		type2.select(GenericInputTypes.GENERIC);
	}
	
	public Input createUIWithData(SADLDescription.Input input) {
		createUI();
		fillWithData(input);
		return this;
	}
	
	private void fillWithData(SADLDescription.Input input) {
		type2.select(input.getType());
		name.setValue(input.getName().getDisplayName());
		description.setValue(getValue(input.getComment()));
		cbMeta.setValue(input.isMeta());
		if (input.isOptional()) {
			optional.select(OPTIONAL);
		} else {
			optional.select(NOT_OPTIONAL);
		}
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
	
	public SADLDescription.Input getSadlInput() {
		SADLDescription.Input input = new SADLDescription.Input();
		input.setName(getNameFromUI(type.getValue().toString()));
		input.setOptional(optional.getValue().equals(OPTIONAL) ? true : false);
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

}
