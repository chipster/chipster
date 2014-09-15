package fi.csc.microarray.client.dataimport;

import java.io.File;

import fi.csc.microarray.databeans.ContentType;

/**
 * Class to represent a single item (file) which will be imported. All 
 * necessary information of the file (input file, output file, action, content 
 * type) are gathered together to this class. Import items are stored to 
 * list in the ImportSession.
 * 
 * 
 * @author mkoski
 *
 */
public class ImportItem {

	/**
	 * An enumeration to describe the action which will be done to file
	 * 
	 * @author mkoski
	 *
	 */
	public enum Action {
		DIRECT("Import directly"), 
		CUSTOM("Use Import tool"), 
		IGNORE("Don't import");
		
		private String name;
		
		private Action(String name){
			this.name = name;
		}
		
		@Override
		public String toString(){
			return name;
		}
	}
	
	private Object input;
	private String inputFilename;
	private ContentType type;
	private Action action;
	/**
	 * Creates a new import item from file. By default the output filename 
	 * is same as input filename.
	 * 
	 * @param input Input file
	 */
	public ImportItem(Object input, String inputFilename, ContentType type) {
		this.input = input;
		this.inputFilename = inputFilename;
		this.type = type;
		this.action = Action.DIRECT; // by default import directly
	}
	
	public Object getInput(){
		return input;
	}

	public Action getAction() {
		return action;
	}
	public String getInputFilename() {
		return inputFilename;
	}

	public ContentType getType() {
		return type;
	}
	public void setAction(Action action){
		this.action = action;
	}

	public void setInputFilename(String inputFilename) {
		this.inputFilename = inputFilename;
	}

	public void setInput(File input) {
		this.input = input;
	}
}
