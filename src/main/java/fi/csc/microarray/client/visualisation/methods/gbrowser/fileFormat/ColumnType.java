package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

/**
 * Possible column types.
 * 
 * @author Petri Klemelä
 * 
 */
public enum ColumnType {
	ID, VALUE, STRAND, SEQUENCE, SKIP, BP_START, BP_END, CHROMOSOME, DESCRIPTION, QUALITY, 
	PARENT_BP_START, PARENT_BP_END, PARENT_ID, PARENT_PART, METADATA, ALLELE, 
	CONSEQUENCE_TO_TRANSCRIPT, POSITION, CIGAR, THICK_START, THICK_END, ITEM_RGB, BLOCK_COUNT, 
	BLOCK_SIZES, BLOCK_STARTS
};