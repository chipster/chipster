package fi.csc.microarray.databeans;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import javax.swing.Icon;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProvider;
import fi.csc.microarray.databeans.features.Modifier;
import fi.csc.microarray.exception.MicroarrayException;

public interface DataManager {

	/**
	 * The initial name for the root folder.
	 */
	public final static String ROOT_NAME = "Datasets";
	
	/**
	 * Adds a listener listening to changes in beans and folders of this manager.
	 */
	public void addDataChangeListener(DataChangeListener listener);
	
	/**
	 * @param enabled if property change events should be sent
	 * @see #addDataChangeListener(DataChangeListener)
	 */
	public void setEventsEnabled(boolean enabled);		

	/**
	 * Returns the root folder, acting as a gateway into the actual data
	 * content under this manager.
	 */
	public DataFolder getRootFolder();
	
	/**
	 * Closes manager and related network connections etc.
	 * Actual implementation depends on the type of manager.
	 */
	public void close();
	
	/**
	 * Creates a folder under this manager. Folder will be created without parent.
	 * 
	 * @param name name for the new folder
	 */
	public DataFolder createFolder(String name);

	/**
	 * Creates a folder under this manager. 
	 * 
	 * @param parent under which folder the new folder is to be created
	 * @param name name for the new folder
	 */
	public DataFolder createFolder(DataFolder root, String name);

	/**
	 * Creates a bean without content, without a parent folder and without sources. If a reference to this bean
	 * is lost it can not be accessed any more.
	 */
	public DataBean createDataBean(String name) throws MicroarrayException;

	/**
	 * Creates a bean with content, without a parent folder and without sources. If a reference to this bean
	 * is lost it can not be accessed any more.
	 */
	public DataBean createDataBean(String name, InputStream content) throws MicroarrayException;

	/**
	 * Creates a bean without content, with a parent folder and with sources.
	 */
	public DataBean createDataBean(String name, DataFolder parent, DataBean... sources) throws MicroarrayException;

	/**
	 * Creates a bean with content, with a parent folder and with sources.
	 */
	public DataBean createDataBean(String name, InputStream content, DataFolder parent, DataBean... sources) throws MicroarrayException;
	
	/**
	 * Guesses MIME content type from a dataset name.
	 */
	public ContentType guessContentType(String name);
	
	/**
	 * Guesses MIME content type from a filename and possibly file content.
	 */
	public ContentType guessContentType(File file);
	
	/**
	 * @return MIME content type for a given extension
	 */
	public ContentType getContentType(String typeName);
	
	/**
	 * Plugs a feature factory, so that it can be used in all beans under this manager.
	 */
	public void plugFeatureFactory(String name, FeatureProvider plugin);
	
	/**
	 * Plugs a modifier (part of Feature API), so that it can be used in all beans under this manager.
	 */
	public void plugModifier(String name, Modifier modifier); 
	public Feature fetchFeature(String featureName, DataBean bean);
	public Modifier fetchModifier(String modifierName);
	
	/**
	 * Plugs a MIME content type, so that it can be used in all beans under this manager.
	 * 
	 * @param mimeType MIME name
	 * @param supported is this a known (supported directly) content type? 
	 * @param description a short textual description
	 * @param extensions file extensions belonging to this type
	 */
	public void plugContentType(String mimeType, boolean supported, boolean binary, String description, Icon icon, String... extensions);
	
	/**
	 * Saves session (all data: beans, folder structure, operation metadata, links etc.) to a file.
	 * File is a zip file with all the data files and one metadata file.
	 * @return count of stored files
	 */
	public void saveSnapshot(File sessionFile, ClientApplication application) throws IOException;

	/**
	 * Load session from a file.
	 * 
	 * @see #saveSnapshot(File, ClientApplication)
	 */
	public List<DataItem> loadSnapshot(File sessionFile, DataFolder parentFolder, ClientApplication application) throws IOException, MicroarrayException;

	
	/**
	 * Load session from an old style session directory (Chipster 1.1 workspace).
	 * 
	 * @see #saveSnapshot(File, ClientApplication)
	 */
	public List<DataItem> loadOldSnapshot(File snapshotDir, DataFolder parentFolder, ClientApplication application) throws IOException, MicroarrayException;

	
	/**
	 * Find and return the first DataItem with the given name.
	 * @param name the name of the DataItem being search for
	 * @return the first found DataItem with given name
	 */
	public DataItem findDataItem(String name);

	/**
	 * Delete DataItem and its children (if any). Root folder cannot be removed.
	 * 
	 * @param data item to be deleted
	 */
	public void delete(DataItem data);

	/**
	 * Return all DataBeans under this manager.
	 */
	public List<DataBean> databeans();

	/**
	 * Return all DataFolders under this manager.
	 */
	public List<DataFolder> folders();

}
