package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

/**
 * FILE_INDEX must not be requested, it's included in results always
 * 
 * @author klemela
 *
 */
public enum ColumnType { ID, VALUE, STRAND, SEQUENCE, SKIP, BP_START, BP_END, CHROMOSOME, 
	DESCRIPTION, FILE_INDEX, QUALITY, PARENT_BP_START, PARENT_BP_END, PARENT_ID, PARENT_PART};