package fi.csc.microarray.client.dataimport;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataimport.ImportItem.Action;

/**
 * A class to store the information of a single import process. The necessary 
 * information are for example files to be imported and destination folder.
 * 
 * @author Mikko Koski
 *
 */
public class ImportSession{
	
	/**
	 * An enumeration to describe the source of imported data in this session.
	 * 
	 * @author Mikko Koski
	 *
	 */
	public enum Source {FILES, FOLDER, URL, CLIPBOARD, TEXT}
	
	private Source source;
	private List<ImportItem> items;
	private String destinationFolder;
	private boolean useSameDescriptions;
	private boolean skipActionChooser;
	
	public ImportSession(Source source, File[] files, String destinationFolder, boolean skipActionChooser) {
		this(source, arrayToList(files), destinationFolder, skipActionChooser);
	}
	
	public ImportSession(Source source, List<File> files, String destinationFolder, boolean skipActionChooser) {
		this.source = source;
		this.items = initializeItems(files);
		this.destinationFolder = destinationFolder;
		this.useSameDescriptions = false;
		this.skipActionChooser = skipActionChooser;
	}

	private static List<ImportItem> initializeItems(List<File> files) {
		List<ImportItem> returnItems = new ArrayList<ImportItem>();
		
		for(int i = 0; i < files.size(); i++){
			ImportItem item = new ImportItem(files.get(i));
			
			// Set default action
			if(ImportUtils.isFileSupported(files.get(i))){
				item.setAction(ImportItem.Action.DIRECT);
			} else {
				item.setAction(ImportItem.Action.CUSTOM);
			}
			
			// Set content type
			item.setType(Session.getSession().getDataManager().guessContentType(files.get(i)));
		
			returnItems.add(item);
		}
		
		return returnItems;
	}

	public List<ImportItem> getImportItems() {
		return items;
	}

	public Source getSource() {
		return source;
	}

	public void setSource(Source source) {
		this.source = source;
	}
	
	public String getDestinationFolder(){
		return this.destinationFolder;
	}
	
	public void setUseSameDescriptions(boolean useSameDescriptions){
		this.useSameDescriptions = useSameDescriptions;
	}
	
	public boolean getUseSameDescriptions(){
		return this.useSameDescriptions;
	}

	public int getItemCount() {
		return this.items.size();
	}
	
	public ImportItem getItemAtIndex(int index){
		if(index >= getItemCount() || index < 0){
			return null;
		} else {
			return items.get(index);
		}
	}
	
	public boolean hasCustomFiles(){
		for(ImportItem item : items){
			if(item.getAction() == ImportItem.Action.CUSTOM){
				return true;
			}
		}
		return false;
	}
	
	public List<File> getInputFiles(){
		ArrayList<File> files = new ArrayList<File>();
		for(ImportItem item : items){
			files.add(item.getInput());
		}
		return files;
	}
	
	public List<File> getCustomFiles(){
		ArrayList<File> files = new ArrayList<File>();
		for(ImportItem item : items){
			if(item.getAction() == Action.CUSTOM){
				files.add(item.getInput());
			}
		}
		return files;
	}

	public void setActionToItemIndex(int index, ImportItem.Action action) {
		this.items.get(index).setAction(action);
	}
	
	public static List<File> arrayToList(File[] files){
		List<File> fileList = new ArrayList<File>();
		for(File file : files){
			fileList.add(file);
		}
		return fileList;
	}

	public boolean isSkipActionChooser() {
		return skipActionChooser;
	}
}
