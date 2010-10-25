package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

/**
 * Column in tab-separated value file. 
 *
 */
public class ColumnDefinition {

	public static final int TAB_DELIM = -1;

	public ColumnType content;
	public Type type;
	public int length; // in bytes
	public long offset; // sum of preceding columns, initialised in Translator constructor

	public ColumnDefinition(ColumnType content, Type type) {
		this(content, type, TAB_DELIM);
	}

	public ColumnDefinition(ColumnType content, Type type, int length) {
		this.content = content;
		this.type = type;
		this.length = length;
	}

}