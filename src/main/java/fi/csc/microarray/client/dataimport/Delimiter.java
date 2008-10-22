package fi.csc.microarray.client.dataimport;

/**
 * Class which represents a delimiter
 * 
 * @author mkoski
 *
 */
public class Delimiter {
	
	/**
	 * Predefined delimiters
	 */
	public static final Delimiter TAB = new Delimiter("Tab", "\t");
	public static final Delimiter SPACE = new Delimiter("Space", " ");
	public static final Delimiter COMMA = new Delimiter("Comma", ",");
	public static final Delimiter SEMICOLON = new Delimiter("Semicolon", ";");
	
	/**
	 * Array of all predefined delimeters
	 */
	private static final Delimiter[] delimiters = {TAB, SPACE, COMMA, SEMICOLON};
	
	/**
	 * Custom delimiter
	 */
	public static final Delimiter CUSTOM = new Delimiter("Custom", null);
	
	private String name;
	private String string;
	
	/**
	 * Constructor. Note that the constructor is private because there 
	 * is static variables to represent the predefined delimiters and 
	 * <code>stringToDelim</code> method to create a new delimiter 
	 * from a string.
	 * 
	 * @param name delimiter name
	 * @param string delimiter string
	 */
	private Delimiter(String name, String string) {
		this.name = name;
		this.string = string;
	}
	
	public String getName(){
		return this.name;
	}
	
	public String toString(){
		return this.string;
	}
	
	public void setDelimiterString(String string){
		this.string = string;
	}
	
	/**
	 * Returns the predefined delimiters
	 * 
	 * @return array of predefined delimeters
	 */
	public static Delimiter[] values(){
		return delimiters;
	}
	
	/**
	 * Creates a delimiter based on given string. If the string is 
	 * same as one in the predefined delimiters, the predefined one is 
	 * returned. Otherwise the method returns a new custom delimiter
	 * 
	 * @param delimString
	 * @return
	 */
	public static Delimiter stringToDelim(String delimString){
		for(Delimiter delim : Delimiter.values()){
			if(delim.toString().equals(delimString)){
				return delim;
			}
		}
		return createCustomDelimiter(delimString);
	}
	
	/**
	 * Creates custom delimiter
	 * 
	 * @param delimiter
	 * @return
	 */
	private static Delimiter createCustomDelimiter(String delimiter){
		return new Delimiter("Custom", delimiter);
	}
	
	/**
	 * Returns if the delimiter is a custom delimiter
	 * 
	 * @param delimString
	 * @return
	 */
	public static boolean isCustom(String delimString){
		for(Delimiter delim : Delimiter.values()){
			if(delim.toString().equals(delimString)){
				return false;
			}
		}
		return true;
	}
}
