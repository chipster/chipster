package fi.csc.microarray.client.dataimport;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataimport.ImportItem.Action;
import fi.csc.microarray.databeans.ContentType;
import fi.csc.microarray.util.IOUtils;

/**
 * A class to store the information of a single import process. The necessary
 * information are for example files to be imported and destination folder.
 * 
 * @author Mikko Koski
 * 
 */
public class ImportSession {

	/**
	 * An enumeration to describe the source of imported data in this session.
	 * 
	 */
	public enum Source {
		FILE, URL, CLIPBOARD
	}

	private Source source;
	private List<ImportItem> items;
	private String destinationFolder;
	private boolean useSameDescriptions;
	private boolean skipActionChooser;

	public ImportSession(Source source, Object[] inputs,
			String destinationFolder, boolean skipActionChooser) {
		
		this(source, Arrays.asList(inputs), destinationFolder, skipActionChooser);
	}

		
	public ImportSession(Source source, List<Object> inputs,
			String destinationFolder, boolean skipActionChooser) {
		this.source = source;
		this.items = toImportItems(inputs);
		this.destinationFolder = destinationFolder;
		this.useSameDescriptions = false;
		this.skipActionChooser = skipActionChooser;
	}

	private List<ImportItem> toImportItems(List<Object> inputs) {
		List<ImportItem> items = new LinkedList<ImportItem>();
		
		for (Object input : inputs) {
			String name = IOUtils.getFilename(input);
			ContentType type = Session.getSession().getDataManager().guessContentType(name);
			items.add(new ImportItem(input, name, type));			
		}
		
		return items;
	}
	
	public void makeLocal(Runnable whenReady) throws IOException {

		if (!isLocal()) {
			for (ImportItem item : items) {
				Object input = item.getInput();
				if (!(input instanceof File)) {
					if (input instanceof URL) {
						URL url = (URL)input;
						File file = ImportUtils.createTempFile(ImportUtils.URLToFilename(url), ImportUtils.getExtension(ImportUtils.URLToFilename(url)));
						ImportUtils.getURLFileLoader().loadFileFromURL(url, file, destinationFolder, whenReady);
						item.setInput(file);
						
					} else {
						throw new RuntimeException("unknown input type: " + input.getClass().getSimpleName());
					}
				}
			}
			source = Source.FILE;
		} else {		
			SwingUtilities.invokeLater(whenReady);
		}
	}
	
	public List<ImportItem> getImportItems() {
		return items;
	}

	public Source getSource() {
		return source;
	}

	public String getDestinationFolder() {
		return this.destinationFolder;
	}

	public void setUseSameDescriptions(boolean useSameDescriptions) {
		this.useSameDescriptions = useSameDescriptions;
	}

	public boolean getUseSameDescriptions() {
		return this.useSameDescriptions;
	}

	public int getItemCount() {
		return this.items.size();
	}

	public ImportItem getItemAtIndex(int index) {
		if (index >= getItemCount() || index < 0) {
			return null;
		} else {
			return items.get(index);
		}
	}

	public boolean hasCustomFiles() {
		for (ImportItem item : items) {
			if (item.getAction() == ImportItem.Action.CUSTOM) {
				return true;
			}
		}
		return false;
	}

	public List<ImportItem> getInputFiles() {
		return items;
	}

	public List<ImportItem> getCustomFiles() {
		ArrayList<ImportItem> files = new ArrayList<ImportItem>();
		for (ImportItem item : items) {
			if (item.getAction() == Action.CUSTOM) {
				files.add(item);
			}
		}
		return files;
	}

	public boolean isSkipActionChooser() {
		return skipActionChooser;
	}


	public boolean isLocal() {
		return this.source == Source.FILE;
	}
}
