package fi.csc.microarray.client.dataimport;

import java.awt.Color;

/**
 * Enumeration of the column types
 * 
 * @author mkoski
 *
 */
public enum ColumnType {

	ROW_NUMBER("Row number", null, null, false),
	IDENTIFIER_LABEL("Identifier", "identifier", new Color(245,242,224), false),
	SAMPLE_LABEL("Sample", "sample", new Color(175,208,175), true),
	SAMPLE_BG_LABEL("Sample BG", "samplebg",  new Color(216,225,202), true),
	CONTROL_LABEL("Control", "control",  new Color(161,182,208), true),
	CONTROL_BG_LABEL("Control BG", "controlbg",  new Color(211,216,213), true),
	FLAG_LABEL("Flag", "flag", new Color(243,232,163), false),
	ANNOTATION_LABEL("Annotation", "annotation", new Color(221,217,202), false),
	UNUSED_LABEL("Unused", null, Color.WHITE, false);

	/**
	 * Column name
	 */
	private String name;
	
	/**
	 * Column identifiers. The R-scripts uses these once.
	 */
	private String title;
	
	/**
	 * Column color when selected on the import preview table
	 */
	private Color color;

	/**
	 * Are values in this column always numeric?
	 */
	private boolean numeric;
	
	private ColumnType(String name, String identifier, Color color, boolean numeric) {
		this.name = name;
		this.title = identifier;
		this.color = color;
		this.numeric = numeric;
	}
	
	@Override
	public String toString() {
		return name;
	}
	
	/**
	 * Gets column identifier. The R-scripts uses this identifier
	 * @return column indentifier
	 */
	public String getTitle() {
		return title;
	}
	
	/**
	 * Column color
	 * @return column color
	 */
	public Color getColor() {
		return color;
	}
	
	public boolean isNumeric() {
		return numeric;
	}
}
