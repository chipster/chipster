package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Lists all supported data types.
 * 
 * @author Petri Klemelä
 * 
 */
public enum DataType {
	ID, 
	VALUE, 
	STRAND, 
	SEQUENCE, 
	START, 
	END, 
	CHROMOSOME, 
	QUALITY, 
	CIGAR, 
	THICK_START, 
	THICK_END, 
	ITEM_RGB, 
	BLOCK_COUNT, 
	BLOCK_SIZES, 
	BLOCK_STARTS, 
	MATE_POSITION, 
	COVERAGE, 
	COVERAGE_ESTIMATE_FORWARD, 
	COVERAGE_ESTIMATE_REVERSE, 
	LOSS, 
	GAIN, 
	FLOAT_LIST, 
	COVERAGE_AVERAGE, 
	CANCEL, 
	REGION 
	};