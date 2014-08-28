package fi.csc.microarray.databeans;

import javax.swing.Icon;

import fi.csc.microarray.constants.VisualConstants;

/**
 * 
 * Content type for DataBeans.
 * 
 * Type field contains the MIME type for the bean.
 * 
 * Does not override equals() or hashCode()!
 * 
 * @author hupponen
 *
 */
public class ContentType {

	
	private String type;
	private boolean supported;
	private String description;
	private String[] extensions;
	private String iconPath;
	private Icon icon;
	private boolean binary;
	
	public ContentType(String name, boolean supported, boolean binary, String description, String iconPath, String... extensions) {
		this.type = name;
		this.supported = supported;
		this.description = description;
		this.iconPath = iconPath;
		this.extensions = extensions;
		this.binary = binary;
	}

	public String getDescription() {
		return description;
	}

	public String[] getExtensions() {
		return extensions;
	}

	public String getType() {
		return type;
	}

	public boolean isSupported() {
		return supported;
	}

	public Icon getIcon() {
		if (icon == null) {
			icon = VisualConstants.getIcon(iconPath);			 
		}
		
		return icon;
	}
	
	public boolean isBinary() {
		return binary;
	}

}

