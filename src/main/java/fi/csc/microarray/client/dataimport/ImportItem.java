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
	
	private File input;
	private File output;
	private Action action;
	private ContentType type;
	
	/**
	 * Creates a new import item from file. By default the output filename 
	 * is same as input filename.
	 * 
	 * @param input Input file
	 */
	public ImportItem(File input) {
		this.input = input;
		
		// Default output name is same as input name
		this.output = input.getAbsoluteFile();
	}
	
	public File getInput(){
		return input;
	}

	public File getOutput() {
		return output;
	}

	public Action getAction() {
		return action;
	}

	public void setFilename(String output) {
		this.output = new File(this.output.getParentFile().getPath() + File.separator + output);
	}

	public ContentType getType() {
		return type;
	}

	public void setType(ContentType type) {
		this.type = type;
	}
	
	public void setAction(Action action){
		this.action = action;
	}
	
	public String toString(){
		return "ImportItem [input: " + input.getName() + ", output: " + output.getName() +
				", type: " + type.getType() + ", action: " + action.toString();
	}
}
