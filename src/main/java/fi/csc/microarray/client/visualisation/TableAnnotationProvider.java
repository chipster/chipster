package fi.csc.microarray.client.visualisation;

import java.util.HashMap;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;

/**
 * Reads table formatted DataBean and looks for common row annotation and description columns. 
 * Found information can be queried using the row identifier (i.e. anonymous first column).
 * If annotations are not available raw identifiers are returned.
 *  
 * @author Aleksi Kallio
 *
 */
public class TableAnnotationProvider {

	private static final String TABLE_COLUMN_QUERY = "/column/";
	private static final String COLUMN_ALL = "*";
	private static final String COLUMN_DESCRIPTION = "description";
	private static final String COLUMN_SYMBOL = "symbol";

	private HashMap<String, String> annotatedIdentifiers = new HashMap<String, String>();
	private HashMap<String, String> descriptions = new HashMap<String, String>();
	
	private DataBean data; //Just to check validity of this instance 
	
	public TableAnnotationProvider(DataBean bean) throws MicroarrayException {
		this(bean.queryFeatures(TABLE_COLUMN_QUERY + COLUMN_ALL).asTable());		
	}

	public TableAnnotationProvider(Table table) throws MicroarrayException {
		while (table.nextRow()) {
			String identifier = table.getStringValue(" ");
			String symbol = table.getStringValue(COLUMN_SYMBOL);
			String actualDescription = table.getStringValue(COLUMN_DESCRIPTION);
			String annotatedIdentifier = symbol != null ? symbol + " (" + identifier + ")" : identifier;
			String description = actualDescription != null ? actualDescription : annotatedIdentifier;
			annotatedIdentifiers.put(identifier, annotatedIdentifier);
			descriptions.put(identifier, description);
		}
	}

	/**
	 * Returns annotated row identifier (i.e. human understandable name for the data point). 
	 * If not available then raw identifier is returned.
	 */
	public String getAnnotatedRowname(String rowIdentifier) {
		return annotatedIdentifiers.get(rowIdentifier);
	}

	/**
	 * Returns row description. If not available then annotation or raw identifier is returned.
	 */
	public String getRowDescription(String rowIdentifier) {
		return descriptions.get(rowIdentifier);
	}
	
	public boolean isForData(DataBean data){
		return this.data != null ? this.data.equals(data) : false;
	}
}
