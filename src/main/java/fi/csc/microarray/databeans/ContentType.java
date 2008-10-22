package fi.csc.microarray.databeans;

import javax.swing.Icon;

/**
 * 
 * Content type for DataBeans.
 * 
 * Type field contains the MIME type for the bean.
 * 
 * 
 * @author hupponen
 *
 */
public class ContentType {

	
	private String type;
	private boolean supported;
	private String description;
	private String[] extensions;
	private Icon icon;
	private boolean binary;
	
	public ContentType(String name, boolean supported, boolean binary, String description, Icon icon, String... extensions) {
		this.type = name;
		this.supported = supported;
		this.description = description;
		this.icon = icon;
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
		return icon;
	}
	
	public boolean isBinary() {
		return binary;
	}

}

