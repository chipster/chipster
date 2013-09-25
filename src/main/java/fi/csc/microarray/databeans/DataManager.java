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
import org.eclipse.jetty.util.IO;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.session.SessionLoader;
import fi.csc.microarray.client.session.SessionLoader.LoadMethod;
import fi.csc.microarray.client.session.SessionSaver;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.StorageMethod;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProvider;
import fi.csc.microarray.databeans.features.Modifier;
import fi.csc.microarray.databeans.handlers.LocalFileDataBeanHandler;
import fi.csc.microarray.databeans.handlers.ZipDataBeanHandler;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.Strings;

public class DataManager {

	private static final String TEMP_DIR_PREFIX = "chipster";
	private static final int MAX_FILENAME_LENGTH = 256;
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
	private HashMap<String, TypeTag> tagMap = new HashMap<String, TypeTag>();
	
	private LinkedList<DataChangeListener> listeners = new LinkedList<DataChangeListener>();
	
	private boolean eventsEnabled = false;

	private DataFolder rootFolder;	
	private File repositoryRoot;

	private ZipDataBeanHandler zipDataBeanHandler = new ZipDataBeanHandler(this);
	private LocalFileDataBeanHandler localFileDataBeanHandler = new LocalFileDataBeanHandler(this);
	
	public DataManager() throws IOException {
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
	 * Find and return all DataBeans with the given name.
	 * @param name the name of the DataBean being search for
	 * @return A list of found DataBeans with given name. Empty list is returned if none was found. 
	 */
	public LinkedList<DataBean> getDataBeans(String name) {
		
		LinkedList<DataBean> list = new LinkedList<>();
		
		for (DataBean dataBean : databeans()) {
			if (dataBean.getName().equals(name)) {
				list.add(dataBean);
			}
		}
		return list;
	}
	
	
	/**
	 * Create a local temporary file DataBean without content, without a parent folder and without sources. 
	 * If a reference to this bean is lost it can not be accessed any more.
	 */
	public DataBean createDataBean(String name) throws MicroarrayException {
		File contentFile;
		try {
			contentFile = createNewRepositoryFile(name);
		} catch (IOException e) {
			throw new MicroarrayException(e);
		}
		
		return createDataBean(name, StorageMethod.LOCAL_TEMP, null, new DataBean[] {}, contentFile);
	}

	/**
	 * Create a local temporary file DataBean with content, without a parent 
	 * folder and without sources. If a reference to this bean
	 * is lost it can not be accessed any more.
	 */
	public DataBean createDataBean(String name, InputStream content) throws MicroarrayException {
		return createDataBean(name, content, null, new DataBean[] {});
	}

	/**
	 * Create a local file DataBean.
	 * The file is used directly, the contents are not copied anywhere.
	 * 
	 */
	public DataBean createDataBean(String name, File contentFile) throws MicroarrayException {		
		return createDataBean(name, StorageMethod.LOCAL_USER, null, new DataBean[] {}, contentFile);
	}

	/**
	 * For now, only file URLs are supported.
	 * 
	 */
	public DataBean createDataBean(String name, URL url) throws MicroarrayException {
		File contentFile;
		try {
			contentFile = new File(url.toURI());
		} catch (Exception e) {
			throw new IllegalArgumentException("Could not convert " + url + " to a file");
		}
		
		return createDataBean(name, StorageMethod.LOCAL_USER, null, new DataBean[] {}, contentFile);
	}

	/**
	 * Create a zip file DataBean. Bean contents are already in the zipFile and can 
	 * be found using the zipEntryName.
	 * 
	 * @param name
	 * @param zipFile
	 * @param zipEntryName
	 * @return
	 * @throws MicroarrayException
	 */
	public DataBean createDataBean(String name, File zipFile, String zipEntryName) throws MicroarrayException {
		URL url;
		try {
			 url = new URL(zipFile.toURI().toURL(), "#" + zipEntryName);
		} catch (MalformedURLException e) {
			throw new MicroarrayException(e);
		}
		
		DataBean dataBean = new DataBean(name, StorageMethod.LOCAL_SESSION, "", url, guessContentType(name), new Date(), new DataBean[] {}, null, this, zipDataBeanHandler);
		dispatchEventIfVisible(new DataItemCreatedEvent(dataBean));
		return dataBean;
	}

	/**
	 * Create a zip file DataBean. Bean contents are already in the zipFile.
	 * 
	 * @param name
	 * @param url location of the zip file, zip entry name as the fragment
	 * @return
	 * @throws MicroarrayException
	 */
	public DataBean createDataBeanFromZip(String name, URL url) throws MicroarrayException {
		DataBean dataBean = new DataBean(name, StorageMethod.LOCAL_SESSION, "", url, guessContentType(name), new Date(), new DataBean[] {}, null, this, zipDataBeanHandler);
		dispatchEventIfVisible(new DataItemCreatedEvent(dataBean));
		return dataBean;
	}

	
	/**
	 * Create a local temporary file DataBean with content, with a parent folder and with sources.
	 */
	private DataBean createDataBean(String name, InputStream content, DataFolder folder, DataBean... sources) throws MicroarrayException {

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
		DataBean bean = createDataBean(name, StorageMethod.LOCAL_TEMP, folder, sources, contentFile);
		return bean;
	}
	

	/**
	 * The file is used directly, the contents are not copied anywhere.
	 * 
	 */
	private DataBean createDataBean(String name, StorageMethod type, DataFolder folder, DataBean[] sources, File contentFile) throws MicroarrayException {
		URL url;
		try {
			 url = contentFile.toURI().toURL();
		} catch (MalformedURLException e) {
			throw new MicroarrayException(e);
		}
		
		DataBean dataBean = new DataBean(name, type, "", url, guessContentType(name), new Date(), sources, folder, this, localFileDataBeanHandler);
		dispatchEventIfVisible(new DataItemCreatedEvent(dataBean));
		return dataBean;
	}
	
	
	/**
	 * Load session from a file.
	 * 
	 * @see #saveSession(File, ClientApplication)
	 */
	public void loadSession(File sessionFile, LoadMethod loadMethod) {
		SessionLoader sessionLoader;
		try {
			sessionLoader = new SessionLoader(sessionFile, loadMethod, this);
			sessionLoader.loadSession();
		} catch (Exception e) {
			e.printStackTrace();
			Session.getSession().getApplication().showDialog("Opening session failed.", "Unfortunately the session could not be opened properly. Please see the details for more information.", Exceptions.getStackTrace(e), Severity.WARNING, true, DetailsVisibility.DETAILS_HIDDEN, null);
			logger.error("loading session failed", e);
		}
	}

	/**
	 * Saves session (all data: beans, folder structure, operation metadata, links etc.) to a file.
	 * File is a zip file with all the data files and one metadata file.
	 * 
	 * @return true if the session was saved perfectly
	 */
	public boolean saveSession(File sessionFile) {
		SessionSaver sessionSaver = new SessionSaver(sessionFile, this);
		boolean metadataValid = false;
		try {
			// save
			metadataValid = sessionSaver.saveSession();
		} catch (Exception e) {
			// save failed, warn about it
			Session.getSession().getApplication().showDialog("Saving session failed.", "Unfortunately your session could not be saved. Please see the details for more information.\n\nIf you have important unsaved datasets in this session, it might be a good idea to export such datasets using the File -> Export functionality.", Exceptions.getStackTrace(e), Severity.WARNING, true, DetailsVisibility.DETAILS_HIDDEN, null);
			return false;
		}

		// check validation, warn if not valid, return false
		if (!metadataValid) {
			// save was successful but metadata validation failed, warn about it
			String validationDetails = sessionSaver.getValidationErrors();
			Session.getSession().getApplication().showDialog("Problem with saving the session.", "All the datasets were saved successfully, but there were troubles with saving the session information about them. This means that there may be problems when trying to open the saved session file later on.\n\nIf you have important unsaved datasets in this session, it might be a good idea to export such datasets using the File -> Export functionality.", validationDetails, Severity.WARNING, true, DetailsVisibility.DETAILS_HIDDEN, null);
			return false;
		}

		return true;
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
	 * The same as saveLightweightSession(), but before saving the session, makes sure
	 * that all data bean contents have been uploaded to cache. Uploads in necessary.
	 * 
	 * @return true if the session was saved perfectly
	 * @throws Exception 
	 */
	public void saveFeedbackSession(File sessionFile) throws Exception {
	
		// upload beans to cache if necessary
		FileBrokerClient fileBroker = Session.getSession().getServiceAccessor().getFileBrokerClient();
		for (DataBean bean : this.databeans()) {
			try {
				bean.getLock().readLock().lock();

				// bean modified, upload
				if (bean.isContentChanged()) {
					bean.setCacheUrl(fileBroker.addFile(bean.getContentByteStream(), bean.getContentLength(), null)); 
					bean.setContentChanged(false);
				} 

				// bean not modified, check cache, upload if needed
				else if (bean.getCacheUrl() != null && !fileBroker.checkFile(bean.getCacheUrl(), bean.getContentLength())){
					bean.setCacheUrl(fileBroker.addFile(bean.getContentByteStream(), bean.getContentLength(), null));
				}

			} finally {
				bean.getLock().readLock().unlock();
			}
		
		}

		// save lightweight session
		SessionSaver sessionSaver = new SessionSaver(sessionFile, this);
		sessionSaver.saveLightweightSession();
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

	/**
	 * FIXME add locking
	 * 
	 * @param bean
	 * @return
	 * @throws IOException 
	 */
	public OutputStream getContentOutputStreamAndLockDataBean(DataBean bean) throws IOException {

		bean.setContentChanged(true);
		
		// Only local temp beans support output, so convert to local temp bean if needed
		if (!bean.getStorageMethod().equals(StorageMethod.LOCAL_TEMP)) {
			this.convertToLocalTempDataBean(bean);
		}
		
		return bean.getHandler().getOutputStream(bean);
	}

	/**
	 * FIXME locks
	 * 
	 * @param bean
	 * @param out
	 * @throws MicroarrayException
	 * @throws IOException
	 */
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
		// convert non local file beans to local file beans
		if (!(bean.getHandler() instanceof LocalFileDataBeanHandler)) {
			this.convertToLocalTempDataBean(bean);
		}
		
		// get the file
		LocalFileDataBeanHandler handler = (LocalFileDataBeanHandler) bean.getHandler();
		return handler.getFile(bean);
	}
	
	
	private void convertToLocalTempDataBean(DataBean bean) throws IOException {
		// FIXME lock bean
		
		// TODO think about that name
		// copy contents to new file
		File newFile = this.createNewRepositoryFile(bean.getName());
		BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(newFile));
		BufferedInputStream in = new BufferedInputStream(bean.getContentByteStream());
		try {
			IOUtils.copy(in, out);
		} finally {
			IOUtils.closeIfPossible(in);
			IOUtils.closeIfPossible(out);
		}
		// update url, type and handler in the bean
		URL newURL = newFile.toURI().toURL();
		
		bean.setContentUrl(newURL);
		bean.setStorageMethod(StorageMethod.LOCAL_TEMP);
		bean.setHandler(localFileDataBeanHandler);
		bean.setContentChanged(true);
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

	public void plugTypeTag(TypeTag typeTag) {
		this.tagMap.put(typeTag.getName(), typeTag);
	}
	
	public TypeTag getTypeTag(String name) {
		return this.tagMap.get(name);
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
		zipDataBeanHandler.closeZipFiles();
	}
	
	public ZipDataBeanHandler getZipDataBeanHandler() {
		return this.zipDataBeanHandler;
	}
}
