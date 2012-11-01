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
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.jms.JMSException;
import javax.swing.Icon;

import org.apache.log4j.Logger;
import org.eclipse.jetty.util.IO;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.session.SessionLoader;
import fi.csc.microarray.client.session.SessionSaver;
import fi.csc.microarray.databeans.DataBean.ContentLocation;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProvider;
import fi.csc.microarray.databeans.features.Modifier;
import fi.csc.microarray.databeans.handlers.ContentHandler;
import fi.csc.microarray.databeans.handlers.LocalFileContentHandler;
import fi.csc.microarray.databeans.handlers.RemoteContentHandler;
import fi.csc.microarray.databeans.handlers.ZipContentHandler;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.filebroker.NotEnoughDiskSpaceException;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;
import fi.csc.microarray.util.Strings;

public class DataManager {

	public static enum StorageMethod {
		
		LOCAL_USER(true, true),
		LOCAL_TEMP(true, true),
		LOCAL_SESSION(true, false),
		REMOTE_CACHED(false, true),
		REMOTE_STORAGE(false, true);
		
		// Groups that describe how fast different methods are to access.
		// Keep these up-to-date when you add methods!
		public static StorageMethod[] LOCAL_FILE_METHODS = {LOCAL_USER, LOCAL_TEMP};
		public static StorageMethod[] REMOTE_FILE_METHODS = {REMOTE_CACHED, REMOTE_STORAGE};
		public static StorageMethod[] OTHER_SLOW_METHODS = {LOCAL_SESSION};
		
		private boolean isLocal;
		private boolean isRandomAccess;
		
		StorageMethod(boolean isLocal, boolean isRandomAccess) {
			this.isLocal = isLocal;
			this.isRandomAccess = isRandomAccess;
		}

		public boolean isLocal() {
			return isLocal;
		}

		public boolean isRandomAccess() {
			return isRandomAccess;
		}
	}

	private static final String TEMP_DIR_PREFIX = "chipster";
	private static final int MAX_FILENAME_LENGTH = 256;
	private static final Logger logger = Logger.getLogger(DataManager.class);

	/**
	 * Reports session validation related problems.
	 */
	public static class ValidationException extends Exception {

		public ValidationException(String validationDetails) {
			// TODO Auto-generated constructor stub
		}
		
	}

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
	private LinkedList<Module> modules;
	
	private ZipContentHandler zipContentHandler = new ZipContentHandler();
	private LocalFileContentHandler localFileContentHandler = new LocalFileContentHandler();
	private RemoteContentHandler remoteContentHandler = new RemoteContentHandler();
	
	public DataManager() throws Exception {
		rootFolder = createFolder(DataManager.ROOT_NAME);

		// initialize repository 		
		repositoryRoot = createRepository();
	}

	public void setRootFolder(DataFolder folder) {
		this.rootFolder = folder;		
	}

	public File getRepository() {
		return repositoryRoot;
	}
	
	/**
	 * Returns the root folder, acting as a gateway into the actual data
	 * content under this manager.
	 */
	public DataFolder getRootFolder() {
		return rootFolder;
	}
	
	
	public boolean isRootFolder(DataFolder folder) {
		return (rootFolder == folder) && (rootFolder != null);
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
		connectChild(folder, root); // events are dispatched from here
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
		String fileName = beanName.replaceAll("[^\\w\\.\\-_]", "");
		if (fileName.length() < 1) {
			fileName = "data";
		} else if (fileName.length() > MAX_FILENAME_LENGTH) {
			fileName = fileName.substring(0, MAX_FILENAME_LENGTH);
		} 
		
		File file = new File(this.repositoryRoot, fileName);
		
		// if file with the beanName already exists, add running number to the name 
		int indexOfDot = fileName.lastIndexOf(".");
		String newFileName = "";
		for (int i = 1; file.exists() && i < Integer.MAX_VALUE; i++) {
			
			// no dot add to end
			if (indexOfDot < 0 ) {
				newFileName = fileName + "-" + i;
			} 

			// add before last dot
			else {
				newFileName = fileName.substring(0, indexOfDot) + "-" + i + fileName.substring(indexOfDot, fileName.length());
			}
			
			file = new File(this.repositoryRoot, newFileName);
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
		File tempRoot = getTempRoot();
		if (!tempRoot.canWrite()) {
			// give up
			throw new IOException("Could not create repository directory.");
		}
		
		String fileName = TEMP_DIR_PREFIX;
		File repository = new File(tempRoot, fileName);
		
		// if directory with that name already exists, add running number 
		boolean repositoryCreated = false;
		for (int i = 1;  !repositoryCreated && i < 1000; i++) {
			repositoryCreated = repository.mkdir();
			if (!repositoryCreated) {
				repository = new File(tempRoot, fileName + "-" + i);
			}
		}

		if (!repositoryCreated) {
			throw new IOException("Could not create repository directory.");
		}
		
		repository.deleteOnExit();
		return repository;
	}

	private File getTempRoot() {
		File tempDir =  new File(System.getProperty("java.io.tmpdir"));

		// check if temp dir is writeable
		if (!tempDir.canWrite()) {
			// try home dir
			tempDir = new File(System.getProperty("user.home"));
			if (!tempDir.canWrite()) {
				// try current working dir
				tempDir = new File(System.getProperty("user.dir"));
			}
		}
		return tempDir;
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
	 * Find and return the first DataBean with the given name.
	 * @param name the name of the DataBean being search for
	 * @return the first found DataBean with given name
	 */
	public DataBean getDataBean(String name) {
		for (DataBean dataBean : databeans()) {
			if (dataBean.getName().equals(name)) {
				return dataBean;
			}
		}
		return null;
	}
	
	
	/**
	 * Create a local temporary file DataBean without content, without a parent folder and without sources. 
	 * If a reference to this bean is lost it can not be accessed any more.
	 */
	public DataBean createLocalTempDataBean(String name) throws MicroarrayException {
		try {
			File contentFile = createNewRepositoryFile(name);
			DataBean bean = createDataBean(name);
			addUrl(bean, StorageMethod.LOCAL_TEMP, contentFile.toURI().toURL());
			return bean;

		} catch (IOException e) {
			throw new MicroarrayException(e);
		}
	}

	
	/**
	 * Creates new DataBean. Infers content type of the created DataBean from the name.
	 * 
	 * @param name name of the DataBean
	 * @return new DataBean that is not connected to a DataFolder
	 */
	public DataBean createDataBean(String name) throws MicroarrayException {
		DataBean data = new DataBean(name, guessContentType(name), this);
		return data;
	}

	/**
	 * Convenience method for creating a local file DataBean. Initialises the DataBean with local file
	 * location. The file is used directly, the contents are not copied anywhere.
	 * 
	 */
	public DataBean createDataBean(String name, File contentFile) throws MicroarrayException {		
		try {
			DataBean bean = createDataBean(name);
			addUrl(bean, StorageMethod.LOCAL_USER, contentFile.toURI().toURL());
			return bean;
			
		} catch (IOException e) {
			throw new MicroarrayException(e);
		}
	}
	
	/**
	 * Convenience method for creating a local temporary file DataBean with content.
	 * Content stream is read into a temp file and location of the file is stored
	 * to DataBean.
	 */
	public DataBean createDataBean(String name, InputStream content) throws MicroarrayException {

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
		DataBean bean = createDataBean(name);
		try {
			addUrl(bean, StorageMethod.LOCAL_TEMP, contentFile.toURI().toURL());
		} catch (MalformedURLException e) {
			throw new MicroarrayException(e);
		}
		return bean;
	}
	

	/**
	 * Load session from a file.
	 * 
	 * @see #saveSession(File, ClientApplication)
	 */
	public void loadSession(File sessionFile, boolean isDataless) throws Exception {
		SessionLoader sessionLoader = new SessionLoader(sessionFile, isDataless, this);
		sessionLoader.loadSession();
	}

	/**
	 * Saves session (all data: beans, folder structure, operation metadata, links etc.) to a file.
	 * File is a zip file with all the data files and one metadata file.
	 * 
	 * @return true if the session was saved perfectly
	 * @throws Exception 
	 */
	public void saveSession(File sessionFile) throws Exception {

		// save session file
		boolean metadataValid = false;
		SessionSaver sessionSaver = new SessionSaver(sessionFile, this);
		metadataValid = sessionSaver.saveSession();

		// check validation
		if (!metadataValid) {
			// save was successful but metadata validation failed, file might be usable
			String validationDetails = sessionSaver.getValidationErrors();
			throw new ValidationException(validationDetails);
		}
	}

	
	/**
	 * Saves lightweight session (folder structure, operation metadata, links etc.) to a file.
	 * Does not save actual data inside databeans.
	 * 
	 * @return true if the session was saved perfectly
	 * @throws Exception 
	 */
	public void saveLightweightSession(File sessionFile) throws Exception {

		SessionSaver sessionSaver = new SessionSaver(sessionFile, this);
		sessionSaver.saveLightweightSession();
	}

	/**
	 * Returns debug print out of current session state.
	 * 
	 * @return print out of session state
	 */
	public String printSession() {
		StringBuffer buffer = new StringBuffer();
		SessionSaver.dumpSession(rootFolder, buffer);
		return buffer.toString();
	}

	public void saveStorageSession(String name) throws Exception {
		URL sessionUrl = Session.getSession().getServiceAccessor().getFileBrokerClient().saveRemoteSession(name);
		SessionSaver sessionSaver = new SessionSaver(sessionUrl, this);
		sessionSaver.saveStorageSession();
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
	
	/**
	 * Remove all DataBeans and DataFolders, except for the root folder.
	 */
	public void deleteAllDataItems() {
		deleteDataFolder(getRootFolder());
	}
	
	private void deleteDataBean(DataBean bean) {

		// remove from operation history
		for (DataBean source : databeans()) {
			// we must iterate all datas because links cannot be trusted (they might have been removed by user)
			OperationRecord operationRecord = source.getOperationRecord();
			if (operationRecord != null) {
				operationRecord.removeInput(bean);
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
			disconnectChild(bean, folder);
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

	public OutputStream getContentOutputStreamAndLockDataBean(DataBean bean) throws IOException {

		// only local temp beans support output, so convert to local temp bean if needed
		ContentLocation tempLocalLocation = bean.getContentLocation(StorageMethod.LOCAL_TEMP);
		if (tempLocalLocation == null) {
			this.convertToLocalTempDataBean(bean);
			tempLocalLocation = bean.getContentLocation(StorageMethod.LOCAL_TEMP);
		}

		// remove all other locations, as they will become obsolete when OutputStream is written to
		while (bean.getContentLocations().size() > 1) {
			for (ContentLocation location : bean.getContentLocations()) {
				if (location != tempLocalLocation) {
					bean.removeContentLocation(location);
					break; // remove outside of the iterator, cannot continue 
				}
			}
		}
		
		return tempLocalLocation.getHandler().getOutputStream(tempLocalLocation); 
	}

	public void closeContentOutputStreamAndUnlockDataBean(DataBean bean, OutputStream out)
			throws MicroarrayException, IOException {
		try {
			out.close();
		} finally {
//			this.lock.writeLock().unlock();
		}
		ContentChangedEvent cce = new ContentChangedEvent(bean);
		this.dispatchEventIfVisible(cce);
	}

	public File getLocalFile(DataBean bean) throws IOException {
		
		ContentLocation location = bean.getContentLocation(StorageMethod.LOCAL_FILE_METHODS);
		
		// convert non local file beans to local file beans
		if (location == null) {
			this.convertToLocalTempDataBean(bean);
			location = bean.getContentLocation(StorageMethod.LOCAL_FILE_METHODS);
		}
		
		// get the file
		LocalFileContentHandler handler = (LocalFileContentHandler) location.getHandler();
		return handler.getFile(location);
	}
	
	
	private void convertToLocalTempDataBean(DataBean bean) throws IOException {

		// copy contents to new file
		File newFile = this.createNewRepositoryFile(bean.getName());
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(newFile));
		BufferedInputStream in = new BufferedInputStream(bean.getContentStream(DataNotAvailableHandling.EXCEPTION_ON_NA));
		try {
			IOUtils.copy(in, out);
		} finally {
			IOUtils.closeIfPossible(in);
			IOUtils.closeIfPossible(out);
		}

		// update url, type and handler in the bean
		URL newURL = newFile.toURI().toURL();
		addUrl(bean, StorageMethod.LOCAL_TEMP, newURL);
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
			disconnectChild(folder, parent);
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

	public Iterable<File> listAllRepositories() {

		LinkedList<File> repositories = new LinkedList<File>();
		
		File tempRoot = getTempRoot();
		
		for (File file: tempRoot.listFiles()) {
			
			if (file.isDirectory() && file.getName().startsWith(TEMP_DIR_PREFIX)) {
				
				String postfix = file.getName().substring(TEMP_DIR_PREFIX.length());
				if ("".equals(postfix) || Strings.isIntegerNumber(postfix)) {
					
					repositories.add(file);
				}
			}
		}
		
		return repositories;
	}

	public void flushSession() {
		zipContentHandler.closeZipFiles();
	}

	public void setModules(LinkedList<Module> modules) {
		this.modules = modules;
	}
	
	public void connectChild(DataItem child, DataFolder parent) {

		// was it already connected?
		boolean wasConnected = child.getParent() != null;

		// connect to this
		child.setParent(parent);

		// add
		parent.children.add(child);

		// add type tags to data beans
		if (child instanceof DataBean) {
			try {
				addTypeTags((DataBean) child);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}
		}

		// dispatch events if needed
		if (!wasConnected) {
			dispatchEvent(new DataItemCreatedEvent(child));
		}
	}

	public void disconnectChild(DataItem child, DataFolder parent) {
		// remove connections
		child.setParent(null);

		// remove
		parent.children.remove(child);

		// dispatch events
		dispatchEvent(new DataItemRemovedEvent(child));
	}


	public void addTypeTags(DataBean data) throws IOException {

		for (Module module : modules) {
			try {
				module.addTypeTags(data);
				
			} catch (MicroarrayException e) {
				throw new RuntimeException(e);
			}
		}
	}
	
	/**
	 * Returns the handler instance of a given StorageMethod. Handler instances are DataManager specific.
	 */
	private ContentHandler getHandlerFor(StorageMethod method) {
		switch (method) {
		
		case LOCAL_SESSION:
			return zipContentHandler;

		case LOCAL_TEMP:
		case LOCAL_USER:
			return localFileContentHandler;
			
		case REMOTE_CACHED:
		case REMOTE_STORAGE:
			return remoteContentHandler;
			
		default:
			throw new IllegalArgumentException("unrecognised method: " + method);	
		
		}
	}

	public void addUrl(DataBean bean, StorageMethod method, URL url) {
		bean.addContentLocation(new ContentLocation(method, getHandlerFor(method), url));
	}

	public void putToStorage(DataBean dataBean) throws Exception {

		// check if content is still available
		if (dataBean.getContentLocations().size() == 0) {
			return; // no content, nothing to put to storage
		}
		
		// check if already in storage
		ContentLocation storageLocation = dataBean.getContentLocation(StorageMethod.REMOTE_STORAGE); 
		if (storageLocation != null && storageLocation.getHandler().isAccessible(storageLocation)) {
			return;
		}
		
		// move from cache to storage, if in cache
		for (ContentLocation cacheLocation : dataBean.getContentLocations(StorageMethod.REMOTE_CACHED)) {
			if (cacheLocation != null && cacheLocation.getHandler().isAccessible(cacheLocation)) {
				
				// move file in filebroker
				URL storageURL = Session.getSession().getServiceAccessor().getFileBrokerClient().moveFileToStorage(cacheLocation.getUrl(), dataBean.getContentLength());
				dataBean.addContentLocation(new ContentLocation(StorageMethod.REMOTE_STORAGE, getHandlerFor(StorageMethod.REMOTE_STORAGE), storageURL));

				// remove cache location(s), because it is now obsolete  
				dataBean.removeContentLocations(StorageMethod.REMOTE_CACHED);
				return;
			}
		}

		// if not in cache, upload to storage
		ContentLocation closestLocation = dataBean.getClosestContentLocation();
		URL storageURL = Session.getSession().getServiceAccessor().getFileBrokerClient().addFile(FileBrokerArea.STORAGE, closestLocation.getHandler().getInputStream(closestLocation), closestLocation.getHandler().getContentLength(closestLocation), null);
		dataBean.addContentLocation(new ContentLocation(StorageMethod.REMOTE_STORAGE, getHandlerFor(StorageMethod.REMOTE_STORAGE), storageURL));		
	}

	/**
	 * 
	 * @param bean
	 * @param progressListener
	 * @return null if no valid location available
	 * @throws NotEnoughDiskSpaceException
	 * @throws FileBrokerException
	 * @throws JMSException
	 * @throws IOException
	 * @throws Exception
	 */
	public URL getURLForCompAndUploadToCacheIfNeeded(DataBean bean, CopyProgressListener progressListener) throws NotEnoughDiskSpaceException, FileBrokerException, JMSException, IOException, Exception { 

		URL url = null;
		try {
			bean.getLock().readLock().lock();

			// upload only if no valid storage or cached location is found
			for (ContentLocation location : bean.getContentLocations(StorageMethod.REMOTE_CACHED, StorageMethod.REMOTE_STORAGE)) {
				if (location.getHandler().isAccessible(location)) {
					url = location.getUrl();
					break;
				}
			}
			// need to upload
			if (url == null) {
				url = Session.getSession().getServiceAccessor().getFileBrokerClient().addFile(FileBrokerArea.CACHE, bean.getContentStream(DataNotAvailableHandling.EXCEPTION_ON_NA), bean.getContentLength(), progressListener);
				bean.removeContentLocations(StorageMethod.REMOTE_CACHED);
				addUrl(bean, StorageMethod.REMOTE_CACHED, url);
			}

		} finally {
			bean.getLock().readLock().unlock();
		}
		return url;
	}

}
