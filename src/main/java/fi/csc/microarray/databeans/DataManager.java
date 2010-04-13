package fi.csc.microarray.databeans;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Date;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.swing.Icon;

import org.apache.log4j.Logger;
import org.mortbay.util.IO;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProvider;
import fi.csc.microarray.databeans.features.Modifier;
import fi.csc.microarray.databeans.handlers.DataBeanHandler;
import fi.csc.microarray.databeans.handlers.LocalFileDataBeanHandler;
import fi.csc.microarray.databeans.handlers.ZipDataBeanHandler;
import fi.csc.microarray.databeans.sessions.OldFSSnapshottingSession;
import fi.csc.microarray.databeans.sessions.SnapshottingSession;
import fi.csc.microarray.exception.MicroarrayException;

public class DataManager {

	private static final Logger logger = Logger.getLogger(DataManager.class);

	
	/**
	 * The initial name for the root folder.
	 */
	public final static String ROOT_NAME = "Datasets";

	
	private Map<String, FeatureProvider> factories = new HashMap<String, FeatureProvider>();
	private Map<String, Modifier> modifiers = new HashMap<String, Modifier>();
	
	/** MIME types for the DataBeans */
	private Map<String, ContentType> contentTypes = new HashMap<String, ContentType>();
	
	/** Mapping file extensions to content types */
	private Map<String, String> extensionMap = new HashMap<String, String>();
	
	private LinkedList<DataChangeListener> listeners = new LinkedList<DataChangeListener>();
	
	private boolean eventsEnabled = false;

	
	
	private DataFolder rootFolder;	
	private File repositoryRoot;
	
	public DataManager() throws IOException {
		rootFolder = createFolder(DataManager.ROOT_NAME);

		// initialize repository 		
		repositoryRoot = createRepository();
	}

	public void setRootFolder(DataFolder folder) {
		this.rootFolder = folder;		
	}

	/**
	 * Returns the root folder, acting as a gateway into the actual data
	 * content under this manager.
	 */
	public DataFolder getRootFolder() {
		return rootFolder;
	}

	/**
	 * Creates a folder under this manager. Folder will be created without parent.
	 * 
	 * @param name name for the new folder
	 */
	public DataFolder createFolder(String name) {
		DataFolder folder = new DataFolder(this, name);
		return folder;
	}

	/**
	 * Creates a folder under this manager. 
	 * 
	 * @param parent under which folder the new folder is to be created
	 * @param name name for the new folder
	 */
	public DataFolder createFolder(DataFolder root, String name) {
		DataFolder folder = new DataFolder(this, name);
		root.addChild(folder); // events are dispatched from here
		return folder;
	}

	/**
	 * Adds a listener listening to changes in beans and folders of this manager.
	 */

	public void addDataChangeListener(DataChangeListener listener) {
		logger.debug("adding DataChangeListener: " + listener);
		if (listener == null) {
			throw new IllegalArgumentException("listener cannot be null");
		}
		listeners.add(listener);
	}

	
	
	
	/**
	 * Creates a new empty file in the repository managed by this FSDataManager.
	 * All the files in this repository should be created by this method.
	 * 
	 * The actual contents of the files may be added either by the 
	 * createDataBean(..., InputStream) methods of this manager, or 
	 * externally and then using the createDataBean(... File) methods
	 * to create the DataBean.
	 * 
	 * This is needed to to avoid overwriting data, in the case of
	 * duplicate DataBean names.
	 * 
	 * @author Taavi Hupponen
	 * 
	 * @param beanName
	 * @return
	 * @throws IOException
	 */
	public synchronized File createNewRepositoryFile(String beanName) throws IOException {
		// FIXME check the file name 
		String fileName = beanName.replaceAll("[^\\wÃ¶Ã¤Ã¥ÃÃÃ\\\\.]", "");
		if (fileName.length() < 1) {
			fileName = "data";
		} else if (fileName.length() > 50) {
			fileName = fileName.substring(0, 50);
		} 
		
		File file = new File(this.repositoryRoot, fileName);
		
		// if file with the beanName already exists, add running number to the name 
		for (int i = 1; file.exists() && i < Integer.MAX_VALUE; i++) {
			file = new File(this.repositoryRoot, fileName + "-" + i);
		}
			
		// create the file
		if (!file.createNewFile()) {
			throw new IOException("Could not create file " + fileName);
		}
		
		// return the file
		file.deleteOnExit();
		return file;
	}

	private File createRepository() throws IOException {
		// get temp dir
		File tempDir =  new File(System.getProperty("java.io.tmpdir"));

		// check if temp dir is writeable
		if (!tempDir.canWrite()) {
			// try home dir
			tempDir = new File(System.getProperty("user.home"));
			if (!tempDir.canWrite()) {
				// try current working dir
				tempDir = new File(System.getProperty("user.dir"));
				if (!tempDir.canWrite()) {
					// give up
					throw new IOException("Could not create repository directory.");
				}
			}
			
		}
		
		
		String fileName = "chipster";
		File repository = new File(tempDir, fileName);
		
		// if file with the beanName already exists, add running number 
		boolean repositoryCreated = false;
		for (int i = 1;  !repositoryCreated && i < 1000; i++) {
			repositoryCreated = repository.mkdir();
			if (!repositoryCreated) {
				repository = new File(tempDir, fileName + "-" + i);
			}
		}

		if (!repositoryCreated) {
			throw new IOException("Could not create repository directory.");
		}
		
		repository.deleteOnExit();
		return repository;
		
	}

	
	
	
	
	/**
	 * @param enabled if property change events should be sent
	 * @see #addDataChangeListener(DataChangeListener)
	 */
	public void setEventsEnabled(boolean enabled) {
		this.eventsEnabled = enabled;		
	}
	
	
	public void dispatchEventIfVisible(DataChangeEvent event) {
		if (event.getDataItem().getParent() != null) {
			dispatchEvent(event);
		}
	}
	
	public void dispatchEvent(DataChangeEvent event) {
		if (eventsEnabled) {
			// dispatch events only for connected datas
			for (DataChangeListener listener : listeners) {
				if (listener == null) {
					logger.error("One of the DataChangeListeners listeners was null.");
				} else {
					logger.debug("Notifying DataChangeListener " + listener.toString());
				}
				try {
					listener.dataChanged(event);
					
				} catch (RuntimeException e) {
					// we will not let GUI problems to stop important DataBean manipulation operations
					// and possibly lead to DataBean model corruption
					logger.error("DataChangeEvent dispatch failed", e);
				}
			}
		}
	}

	
	public DataBean createDataBean(String name, DataFolder parent, DataBean source) throws MicroarrayException {
		return createDataBean(name, parent, wrapSource(source));
	}
	
	public static DataBean[] wrapSource(DataBean source) {
		DataBean[] sources = null;
		
		if (source != null) {
			sources = new DataBean[1];
			sources[0] = source;
		} else {
			sources = new DataBean[0];
		}
		
		return sources;
		
	}


	
	





	
	/**
	 * Guess the MIME content type using the filename.
	 * 
	 * For now, simply use the extension to figure out the mime type.
	 * 
	 * Types are plugged at ApplicationConstants.
	 * 
	 */
	public ContentType guessContentType(String name) {
		
		ContentType type = null;
		if(name.lastIndexOf(".") != -1){
			String extension = name.substring(name.lastIndexOf(".") + 1, name.length()).toLowerCase();
			String typeName = extensionMap.get(extension);
			if (typeName != null) {
				type = contentTypes.get(typeName);
			}
		} 
		
		if (type == null) {
			type = contentTypes.get("application/octet-stream");
		}
		return type;
	}
	
	/**
	 * Guesses MIME content type from a filename and possibly file content.
	 */
	public ContentType guessContentType(File file) {
		return guessContentType(file.getName());
	}
	

	/**
	 * @return MIME content type for a given extension
	 */
	public ContentType getContentType(String typeName) {
		return contentTypes.get(typeName);
	}
	
	/**
	 * Plugs a MIME content type, so that it can be used in all beans under this manager.
	 * 
	 * @param mimeType MIME name
	 * @param supported is this a known (supported directly) content type? 
	 * @param description a short textual description
	 * @param extensions file extensions belonging to this type
	 */
	public void plugContentType(String mimeType, boolean supported, boolean binary, String description, Icon icon, String... extensions) {
		// create the content type
		contentTypes.put(mimeType, new ContentType(mimeType, supported, binary, description, icon, extensions));
		
		
		// add extensions to search map
		for (String extension: extensions) {
			extensionMap.put(extension, mimeType);
		}
	}

	/**
	 * Plugs a modifier (part of Feature API), so that it can be used in all beans under this manager.
	 */
	public void plugModifier(String name, Modifier modifier) {
		modifiers.put(name, modifier);
	}

	/**
	 * Plugs a feature factory, so that it can be used in all beans under this manager.
	 */
	public void plugFeatureFactory(String name, FeatureProvider plugin) {
		logger.debug("plugged " + plugin.getClass().getSimpleName() + " at " + name);
		plugin.setName(name);
		factories.put(name, plugin);
	}


	public Modifier fetchModifier(String modifierName) {	
		return modifiers.get(modifierName);
	}
	
	public Feature fetchFeature(String featureName, DataBean bean) {
		String bestMatch = null; 		
		for (String feature : factories.keySet()) {
			if (featureName.startsWith(feature)) {
				if (bestMatch == null || feature.length() > bestMatch.length()) {
					// current best match
					bestMatch = feature;
				}
			}
		}
		FeatureProvider factory = factories.get(bestMatch);
		if (factory == null) {
			throw new RuntimeException("no feature factory plugged in for \"" + featureName + "\" (total of " + factories.size() + " factories plugged)");
		}
		logger.debug("best match for " + featureName + " was " + (factory != null ? factory.getName() : factory));
		String namePostfix = getNamePostfix(featureName, factory.getName());		
		return factory.createFeature(namePostfix, bean);
	}

	
	
	


	

	
	/**
	 * Find and return the first DataItem with the given name.
	 * @param name the name of the DataItem being search for
	 * @return the first found DataItem with given name
	 */
	public DataItem findDataItem(String name) {
		return findDataItem(name, getRootFolder());
	}
	
	private DataItem findDataItem(String name, DataItem root) {
		DataItem matchingItem = null;
		
		// root item matches
		if (root.getName().equals(name)) {
			return root;
		} 
		
		// root is a folder, search children
		else if (root instanceof DataFolder) {
			for (DataItem child: ((DataFolder)root).getChildren()) {
				matchingItem = findDataItem(name, child);
				if (matchingItem != null) {
					return matchingItem;
				}
			}
		} 
		
		// no match found
		return null;
	}

	
	
	



	
	
	/**
	 * Creates a bean without content, without a parent folder and without sources. If a reference to this bean
	 * is lost it can not be accessed any more.
	 */
	public DataBean createDataBean(String name) throws MicroarrayException {
		return createDataBean(name, null, new DataBean[] {});
	}

	/**
	 * Creates a bean without content, with a parent folder and with sources.
	 */
	public DataBean createDataBean(String name, DataFolder folder, DataBean... sources) throws MicroarrayException {
		File contentFile;
		try {
			contentFile = createNewRepositoryFile(name);
		} catch (IOException e) {
			throw new MicroarrayException(e);
		}
		
		return createDataBean(name, folder, sources, contentFile);
	}
	
	/**
	 * Creates a bean with content, without a parent folder and without sources. If a reference to this bean
	 * is lost it can not be accessed any more.
	 */
	private DataBean createDataBean(String name, InputStream content) throws MicroarrayException {
		return createDataBean(name, content, null, new DataBean[] {});
	}
	
	/**
	 * Creates a bean with content, with a parent folder and with sources.
	 */
	public DataBean createDataBean(String name, InputStream content, DataFolder folder, DataBean... sources) throws MicroarrayException {

		// copy the data from the input stream to the file in repository
		File contentFile;
		try {
			contentFile = createNewRepositoryFile(name);
			InputStream input = new BufferedInputStream(content);
			OutputStream output = new BufferedOutputStream(new FileOutputStream(contentFile));
			IO.copy(input, output);
			input.close();
			output.flush();
			output.close();
		} catch (IOException ioe) {
			throw new MicroarrayException(ioe);
		}

		// create and return the bean
		DataBean bean = createDataBean(name, folder, sources, contentFile);
		return bean;
	}

	
	
	public DataBean createDataBean(String name, File zipFile, String zipEntryName) throws MicroarrayException {
		URL url;
		try {
			 url = new URL(zipFile.toURI().toURL(), "#" + zipEntryName);
		} catch (MalformedURLException e) {
			throw new MicroarrayException(e);
		}
		
		DataBeanHandler handler = new ZipDataBeanHandler();
		DataBean dataBean = new DataBean(name, url, guessContentType(name), new Date(), new DataBean[] {}, null, this, handler);
		dispatchEventIfVisible(new DataItemCreatedEvent(dataBean));
		return dataBean;
	}

	
	public DataBean createLocalFileDataBean(String name, InputStream content) throws MicroarrayException {
		return createDataBean(name, content);
	}
	
	
	/**
	 * The file is used directly, the contents are not copied anywhere.
	 * 
	 */
	public DataBean createDataBean(String name, File contentFile) throws MicroarrayException {		
		return createDataBean(name, null, new DataBean[] {}, contentFile);
	}

	/**
	 * The file is used directly, the contents are not copied anywhere.
	 * 
	 */
	private DataBean createDataBean(String name, DataFolder folder, DataBean[] sources, File contentFile) throws MicroarrayException {
		URL url;
		try {
			 url = contentFile.toURI().toURL();
		} catch (MalformedURLException e) {
			throw new MicroarrayException(e);
		}
		
		DataBeanHandler handler = new LocalFileDataBeanHandler();
		DataBean dataBean = new DataBean(name, url, guessContentType(name), new Date(), sources, folder, this, handler);
		dispatchEventIfVisible(new DataItemCreatedEvent(dataBean));
		return dataBean;
	}
	
	

	/**
	 * Load session from a file.
	 * 
	 * @see #saveSnapshot(File, ClientApplication)
	 */
	public List<DataItem> loadSnapshot(File sessionFile, DataFolder parentFolder, ClientApplication application) throws IOException, MicroarrayException {
		SnapshottingSession session = new SnapshottingSession(this, application);
		List<DataItem> newItems = session.loadFromSnapshot(sessionFile, parentFolder);
		return newItems;
	}

	/**
	 * Load session from an old style session directory (Chipster 1.1 workspace).
	 * 
	 * @see #saveSnapshot(File, ClientApplication)
	 */
	public List<DataItem> loadOldSnapshot(File snapshotDir, DataFolder parentFolder, ClientApplication application) throws IOException, MicroarrayException {
		OldFSSnapshottingSession session = new OldFSSnapshottingSession(this);
		List<DataItem> newItems = session.loadFromSnapshot(snapshotDir, parentFolder);
		return newItems;
	}


	/**
	 * Saves session (all data: beans, folder structure, operation metadata, links etc.) to a file.
	 * File is a zip file with all the data files and one metadata file.
	 * @return count of stored files
	 */
	public void saveSnapshot(File snapshotDir, ClientApplication application) throws IOException {
		SnapshottingSession session = new SnapshottingSession(this, application);
		session.saveSnapshot(snapshotDir);
	}

	/**
	 * Delete DataItem and its children (if any). Root folder cannot be removed.
	 * 
	 * @param data item to be deleted
	 */
	public void delete(DataItem data) {
		
		if (data instanceof DataFolder) {
			deleteDataFolder((DataFolder)data);
			
		} else {
			deleteDataBean((DataBean)data);
		}		
	}
	
	private void deleteDataBean(DataBean bean) {

		// remove from operation history
		for (DataBean source : databeans()) {
			// we must iterate all datas because links cannot be trusted (they might have been removed by user)

			boolean isDirty = false;
			List<DataBinding> bindings = source.getOperation().getBindings();

			if (bindings != null) {
				for (DataBinding binding : bindings) {
					if (binding.getData() == bean) {
						// this operation would become dirty after removing the data
						isDirty = true;
						break;
					}
				}
			}

			if (isDirty) {
				source.getOperation().clearBindings();
			}
		}
		
		// remove links
		for (Link linkType : Link.values()) {
			// Remove outgoing links
			for (DataBean target : bean.getLinkTargets(linkType)) {
				bean.removeLink(linkType, target);
			}
			// Remove incoming links
			for (DataBean source : bean.getLinkSources(linkType)) {
				source.removeLink(linkType, bean);
			}
		}

		// remove bean
		DataFolder folder = bean.getParent();
		if (folder != null) {
			folder.removeChild(bean);
		}
		
		// remove physical file
		bean.delete();
	}

	/**
	 * Return all DataBeans under this manager.
	 */
	public List<DataBean> databeans() {
		
		LinkedList<DataBean> databeans = new LinkedList<DataBean>();
		for (DataFolder folder : folders()) {
			for (DataItem child : folder.getChildren()) {
				if (child instanceof DataBean) {
					databeans.add((DataBean) child);
				}
			}
		}
		return databeans;		
	}

	/**
	 * Return all DataFolders under this manager.
	 */
	public List<DataFolder> folders() {
		return folders(getRootFolder());
	}

	public List<DataFolder> folders(DataFolder parent) {
		LinkedList<DataFolder> folders = new LinkedList<DataFolder>();
		folders.add(parent);
		for (DataItem child : parent.getChildren()) {
			if (child instanceof DataFolder) {
				folders.addAll(folders((DataFolder) child));
			}
		}
		return folders;
	}


	private void deleteDataFolder(DataFolder folder) {

		// remove children
		Iterable<DataItem> children = folder.getChildren();

		// make a copy of the children list to avoid concurrent modification
		List<DataItem> childrenToBeRemoved = new LinkedList<DataItem>();
		for (DataItem item : children) {
			childrenToBeRemoved.add(item);
		}

		// remove all children (recursively)
		for (DataItem item : childrenToBeRemoved) {
			delete(item); 
		}

		// remove this folder (unless root)
		DataFolder parent = folder.getParent();
		if (parent != null) {
			parent.removeChild(folder);
		}
	}

	
	
	private String getNamePostfix(String featureName, String factoryName) {
		if (factoryName.length() > featureName.length()) {
			return "";
		} else {
			String npf = featureName.substring(factoryName.length());
			if (npf.startsWith("/")) {
				return npf.substring(1);
			} else {
				return npf;
			}
		}
	}

}
