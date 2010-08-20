package fi.csc.microarray.databeans;

/**
 * Type tags are attached to databeans. Each of them tells about one aspect 
 * of bean's motivation: what it is meant for and how it should be handled. 
 * The complete type of a bean is made up from the MIME content type and 
 * the set of type tags.  
 * 
 * @author Aleksi Kallio
 *
 * @see DataBean
 */
public class TypeTag {

	private String name;

	public TypeTag(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}
	
}
