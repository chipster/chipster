package fi.csc.microarray.client.dataimport;

import fi.csc.microarray.client.dataimport.trimmer.DataTrimmer;

/**
 * The most important methods that import screen should implement
 * 
 * @author mkoski
 *
 */
public interface ImportScreenModel {

	/**
	 * Gets the conversion model
	 * 
	 * @return conversion model
	 */
	public ConversionModel getConversionModel();
	
	/**
	 * Gets the column type manager
	 * 
	 * @return column type manager
	 */
	public ColumnTypeManager getColumnTypeManager();
	
	/**
	 * Gets the data trimmer
	 * 
	 * @return data trimmer
	 */
	public DataTrimmer getDataTrimmer();
	
	/**
	 * Flag trimmer
	 * 
	 * @return flag trimmer
	 */
	public DataTrimmer getFlagTrimmer();
	
//	/**
//	 * Sets the input file
//	 * 
//	 * @param input input file
//	 */
//	public void setInput(File input);
	
	public void setImportSession(ImportSession importSession);
}
