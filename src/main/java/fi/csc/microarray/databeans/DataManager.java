package fi.csc.microarray.databeans;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import javax.jms.JMSException;

import org.apache.log4j.Logger;
import org.eclipse.jetty.util.IO;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProvider;
import fi.csc.microarray.databeans.features.Modifier;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.databeans.handlers.ContentHandler;
import fi.csc.microarray.databeans.handlers.LocalFileContentHandler;
import fi.csc.microarray.databeans.handlers.RemoteContentHandler;
import fi.csc.microarray.databeans.handlers.ZipContentHandler;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.ChecksumException;
import fi.csc.microarray.filebroker.ChecksumInputStream;
import fi.csc.microarray.filebroker.ContentLengthException;
import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.filebroker.NotEnoughDiskSpaceException;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.security.CryptoKey;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;
import fi.csc.microarray.util.Strings;
import fi.csc.microarray.util.ThreadUtils;

public class DataManager {

	/**
	 * Defines a location where the actual file content of the DataBean can be accessed from.
	 * Location might be local, remote or contained within another file (zip). 
	 */
	public static class ContentLocation {
		
		StorageMethod method;
		URL url;
		private ContentHandler handler;
		
		public ContentLocation(StorageMethod method, ContentHandler handler, URL url) {
			this.method = method;
			this.handler = handler;
			this.url = url;
		}
		
		public StorageMethod getMethod() {
			return method;
		}

		public URL getUrl() {
			return url;
		}

		public ContentHandler getHandler() {
			return handler;
		}
	}

	
	public static enum StorageMethod {
		
		LOCAL_ORIGINAL(true, true),
		LOCAL_TEMP(true, true),
		LOCAL_SESSION_ZIP(true, false),
		REMOTE_ORIGINAL(false, true);
		
		// Groups that describe how fast different methods are to access.
		// Keep these up-to-date when you add methods!
		public static StorageMethod[] LOCAL_FILE_METHODS = {LOCAL_ORIGINAL, LOCAL_TEMP};
		public static StorageMethod[] REMOTE_FILE_METHODS = {REMOTE_ORIGINAL};
		public static StorageMethod[] OTHER_SLOW_METHODS = {LOCAL_SESSION_ZIP};
		
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
		
		/**
		 * Returns value of given name, so that old names are first converted to new naming scheme and
		 * then used to call {@link Enum#valueOf(Class, String)}.
		 * 
		 * @param name name or old name of StorageMethod enum value
		 * 
		 * @return StorageMethod enum value
		 */
		public static StorageMethod valueOfConverted(String name) {
			if ("REMOTE_STORAGE".equals(name)) {
				name = "REMOTE_ORIGINAL";
			} else if ("LOCAL_USER".equals(name)) {
				name = "LOCAL_ORIGINAL";
			} else if ("LOCAL_SESSION".equals(name)) {
				name = "LOCAL_SESSION_ZIP";
			}

			return valueOf(name);
		}
		
	}
	
	private static final String AT_LEAST_ROWS_CACHENAME = "at-least-rows";
	public static final int MAX_ROWS_TO_COUNT = 1000; 
	public static final int MAX_BYTES_TO_COUNT = 100*1024;

	public static final String DATA_NA_INFOTEXT = "Data currently not available";
	private static final String TEMP_DIR_PREFIX = "chipster";
	private static final int MAX_FILENAME_LENGTH = 256;
	private static final Logger logger = Logger.getLogger(DataManager.class);

	/**
	 * Reports session validation related problems.
	 */
	@SuppressWarnings("serial")
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
	
	// by default there are max 5 simultaneous http connections, but type tagging
	// creates also JMS traffic
	private ExecutorService executor = Executors.newFixedThreadPool(10);
	
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
		// this could be replaced with Java 7's Files.createTempDirectory
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
	public void plugContentType(String mimeType, boolean supported, boolean binary, String description, String iconPath, String... extensions) {
		// create the content type
		contentTypes.put(mimeType, new ContentType(mimeType, supported, binary, description, iconPath, extensions));
		
		
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
	public DataBean createLocalTempDataBean(String name) throws MicroarrayException {
		try {
			File contentFile = createNewRepositoryFile(name);
			DataBean bean = createDataBean(name);
			addContentLocationForDataBean(bean, StorageMethod.LOCAL_TEMP, contentFile.toURI().toURL());
			return bean;

		} catch (IOException | ContentLengthException e) {
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
	 * Convenience method for creating a filebroker DataBean. Infers content type of the created DataBean from the name and
	 * tries to get the content length of the file from filebroker.
	 * 
	 * @param name name of the DataBean
	 * @param updateContentLength 
	 * @return new DataBean that is not connected to a DataFolder
	 */
	public DataBean createDataBean(String name, String dataId, boolean updateContentLength) throws MicroarrayException {
		DataBean data = new DataBean(name, guessContentType(name), this, dataId);
		if (updateContentLength) {
			/* Update content length. Usually it is updated when content location 
			 * is added, but filebroker is not a content location.  
			 */
			getContentLength(data);
		}
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
			addContentLocationForDataBean(bean, StorageMethod.LOCAL_ORIGINAL, contentFile.toURI().toURL());
			return bean;
			
		} catch (IOException | ContentLengthException e) {
			throw new MicroarrayException(e);
		}
	}
	
	/**
	 * Convenience method for creating a remote url DataBean. Initialises the DataBean with the url
	 * location. The url is used directly, the contents are not copied anywhere.
	 * @throws IOException 
	 * @throws ContentLengthException 
	 * 
	 */
	public DataBean createDataBean(String name, URL url) throws MicroarrayException, ContentLengthException, IOException {		
		DataBean bean = createDataBean(name);
		
		addContentLocationForDataBean(bean, StorageMethod.REMOTE_ORIGINAL, url);
		
		return bean;
	}
	
	/**
	 * Convenience method for creating a local temporary file DataBean with content.
	 * Content stream is read into a temp file and location of the file is stored
	 * to DataBean.
	 * @throws IOException 
	 */
	public DataBean createDataBean(String name, InputStream content) throws MicroarrayException, IOException {

		// copy the data from the input stream to the file in repository
		File contentFile;
		contentFile = createNewRepositoryFile(name);
		InputStream input = new BufferedInputStream(content);
		OutputStream output = new BufferedOutputStream(new FileOutputStream(contentFile));
		try {
			IO.copy(input, output);
			output.flush();
		} finally {
			IOUtils.closeIfPossible(input);
			IOUtils.closeIfPossible(output);
		}

		// create and return the bean
		DataBean bean = createDataBean(name);
		try {			
			addContentLocationForDataBean(bean, StorageMethod.LOCAL_TEMP, contentFile.toURI().toURL());			
		} catch (MalformedURLException e) {
			throw new MicroarrayException(e);
		} catch (ContentLengthException e) {
			// shouldn't happen, because newly created bean doesn't have size set
			logger.error(e,e);
		}
		return bean;
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

	public ChecksumInputStream getContentStream(DataBean bean, DataNotAvailableHandling naHandling) throws IOException {		
		
		/* It's convenient to return a ChecksumInputStream, because it transparently handles situation of disabled 
		 * checksum calculation. Method getContentStreamImpl creates many types of streams, so wrap a ChecksumInputStream
		 * around the returned stream if necessary. 
		 */
		InputStream baseStream  = getBaseContentStream(bean, naHandling);
		
		if (baseStream instanceof ChecksumInputStream) {
			// JMSFileBrokerClient enables checksum calculation by default 			
			return (ChecksumInputStream) baseStream;
		} else if (baseStream == null) {
			return null;
		} else {
			/* Disable checksum calculation for other streams (local copies etc.), 
			 * because of supposed performance implications. Nevertheless, enabling
			 * might work as well and would provide more thorough data integrity check.
			 */
			return new ChecksumInputStream(baseStream, false);
		}
	}
		
	private InputStream getBaseContentStream(DataBean bean, DataNotAvailableHandling naHandling) throws IOException {

		// try local content locations first
		ContentLocation location = getClosestContentLocation(bean);
		
		
		if (location != null) {
			
			// local available TODO maybe check if it really is available
			return location.getHandler().getInputStream(location);
		 
		} 
		
		// try from filebroker
		Exception remoteException;
		try {
			return Session.getSession().getServiceAccessor().getFileBrokerClient().getInputStream(bean.getId());
		} catch (Exception e) {
			remoteException = e;
		}
		
		// not available
		switch (naHandling) {
		case EMPTY_ON_NA:
			return new ByteArrayInputStream(new byte[] {});

		case INFOTEXT_ON_NA:
			return new ByteArrayInputStream(DATA_NA_INFOTEXT.getBytes());

		case NULL_ON_NA:
			return null;

		default:
			String message = "data contents not available";
			if (remoteException != null) {
				message += ": " + remoteException.getMessage(); 
			}
			throw new RuntimeException(message, remoteException);	
		}

	}

	
	/**
	 * A convenience method for gathering streamed binary content into
	 * a byte array. Returns null if none of the content locations are available. 
	 * 
	 * @throws IOException 
	 * 
	 * @see #getContentStream()
	 */
	public byte[] getContentBytes(DataBean bean, DataNotAvailableHandling naHandling) throws IOException {
		return getContentBytes(bean, -1, naHandling); // -1 means "no max length"
	}

	/**
	 * A convenience method for gathering streamed binary content into
	 * a byte array. Gathers only maxLength first bytes. Returns null if 
	 * none of the content locations are available. 
	 * 
	 * @throws IOException 
	 * 
	 * @see #getContentStream()
	 */
	public byte[] getContentBytes(DataBean bean, long maxLength, DataNotAvailableHandling naHandling) throws IOException {
		
		InputStream in = null;
		try {
			in = getContentStream(bean, naHandling);
			if (in != null) {
				return Files.inputStreamToBytes(in, maxLength);
				
			} else {
				return null;
			}
			
		} finally {
			IOUtils.closeIfPossible(in);
		}
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
		
		// change data id
		bean.setId(CryptoKey.generateRandom());
		bean.setChecksum(null);
		bean.setSize(null);
		
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
	
	/**
	 * Get a local file with random access support, in practice a temp file. If random access support 
	 * isn't really needed, please use method getLocalFile().
	 * 
	 * 
	 * @param bean
	 * @return
	 * @throws IOException
	 */
	public File getLocalRandomAccessFile(DataBean bean) throws IOException {
		
		//Check if there is a suitable location already
		ContentLocation location = getClosestRandomAccessContentLocation(bean);		
		if (location == null || !location.getMethod().isLocal()) {
			//Closest random access location isn't local, make a local copy (temp files do support random access)
			this.convertToLocalTempDataBean(bean);
		}
				
		//Now this is a local random access copy
		location = getClosestRandomAccessContentLocation(bean);
		
		// get the file
		LocalFileContentHandler handler = (LocalFileContentHandler) location.getHandler();
		return handler.getFile(location);
	}

	
	/**
	 * <p>Returns content size in bytes. Returns -1 if 
	 * data is not available.</p>
	 * 
	 * <p>Size reported by content location is stored in DataBean in
	 * read from there when possible.</p>
	 */
	public Long getContentLength(DataBean bean) {
		if (bean.getSize() == null) {
			try {
				ContentLocation location = getClosestContentLocation(bean);
				if (location != null) {
					bean.setSize(getContentLength(location));
				} else {
					Long size = Session.getSession().getServiceAccessor().getFileBrokerClient().getContentLength(bean.getId());
					if (size != null) {
						bean.setSize(size);
					} else {
						return -1l;
					}
				}
			} catch (Exception e) {
				throw new RuntimeException(e);
			}
		}
		return bean.getSize();
	}
	
	public Long getContentLength(ContentLocation location) throws IOException {
		return location.getHandler().getContentLength(location);
	}
	
	
	private void convertToLocalTempDataBean(DataBean bean) throws IOException {
		
		try {
			// copy contents to new file
			File newFile = this.createNewRepositoryFile(bean.getName());
			BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(newFile));
			ChecksumInputStream inputStream = getContentStream(bean, DataNotAvailableHandling.EXCEPTION_ON_NA);
			BufferedInputStream in = new BufferedInputStream(inputStream);
			try {
				IOUtils.copy(in, out);

				setOrVerifyChecksum(bean, inputStream.verifyChecksums());
			} finally {
				IOUtils.closeIfPossible(in);
				IOUtils.closeIfPossible(out);
			}

			// update url, type and handler in the bean
			URL newURL = newFile.toURI().toURL();
			addContentLocationForDataBean(bean, StorageMethod.LOCAL_TEMP, newURL);
		} catch (ChecksumException | ContentLengthException e) {
			// corrupted data  
			throw new IOException();
		}
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
	
	public void connectChildren(final List<? extends DataItem> children, final DataFolder parent) {
		
		final ArrayList<AddTypeTagsCallable> callables = new ArrayList<>();
		final ArrayList<DataItemCreatedEvent> events = new ArrayList<>();
		
		// edit databeans in EDT
		ThreadUtils.runInEDT(new Runnable() {					
			@Override
			public void run() {
				for (DataItem child : children) {
					// was it already connected?
					boolean wasConnected = child.getParent() != null;
					if (!wasConnected) {
						// prepare event, but don't send it yet
						events.add(new DataItemCreatedEvent(child));
					}
					
					// connect to this
					child.setParent(parent);

					// add
					parent.children.add(child);

					// prepare type tagging callables
					if (child instanceof DataBean) {
						callables.add(new AddTypeTagsCallable((DataBean) child));
					}
				}

			}
		});
			
		// run type tagging in parallel in background threads
		try {
			// run callables and wait until all have finished
			List<Future<Object>> futures = executor.invokeAll(callables);
			
			for (Future<Object> future : futures) {
				// check callable for exception, and throw if found
				future.get();
			}
			
			// dispatch events in EDT
			ThreadUtils.runInEDT(new Runnable() {					
				@Override
				public void run() {
					for (DataItemCreatedEvent event : events) {
						dispatchEvent(event);
					}
				}
			});
			
		} catch (InterruptedException | ExecutionException e) {
			Session.getSession().getApplication().reportExceptionThreadSafely(e);
		}
	}
	
	public class AddTypeTagsCallable implements Callable<Object> {
		private DataBean child;

		public AddTypeTagsCallable(DataBean child) {
			this.child = child;
		}

		@Override
		public Object call() throws Exception {
			addTypeTagsOfEachModule(child);
			return null;
		}
	}
	
	public void connectChild(DataItem child, DataFolder parent) {
		ArrayList<DataItem> list = new ArrayList<>();
		list.add(child);
		connectChildren(list, parent);
	}

	public void disconnectChild(DataItem child, DataFolder parent) {
		// remove connections
		child.setParent(null);

		// remove
		parent.children.remove(child);

		// dispatch events
		dispatchEvent(new DataItemRemovedEvent(child));
	}


	public void addTypeTagsOfEachModule(DataBean data) throws IOException {

		if (!data.isTagsSet()) {
			for (Module module : modules) {
				try {
					module.addTypeTags(data);

				} catch (MicroarrayException e) {
					throw new RuntimeException(e);
				}
			}
		}
		data.setTagsSet(true);
	}
	
	/**
	 * Returns the handler instance of a given StorageMethod. Handler instances are DataManager specific.
	 */
	private ContentHandler getHandlerFor(StorageMethod method) {
		switch (method) {
		
		case LOCAL_SESSION_ZIP:
			return zipContentHandler;

		case LOCAL_TEMP:
		case LOCAL_ORIGINAL:
			return localFileContentHandler;
					
		case REMOTE_ORIGINAL:
			return remoteContentHandler;
			
		default:
			throw new IllegalArgumentException("unrecognised method: " + method);	
		
		}
	}
	
	public void addContentLocationForDataBean(DataBean bean, StorageMethod method, URL url) throws ContentLengthException, IOException {
		ContentLocation location = new ContentLocation(method, getHandlerFor(method), url);
		try {
			setOrVerifyContentLength(bean, getContentLength(location));
		} catch (IOException e) {
			throw new IOException("content length not available: " + e);
		} catch (ContentLengthException e) {
			try {
				DataManager manager = Session.getSession().getDataManager();
				String msg;
				msg = "Wrong content length for dataset " + bean.getName() + ". "
						+ " In ContentLocation " + location.getUrl() +  ", length is " + getContentLength(location) + " bytes. ";
				msg += "Content locations: ";
				for (ContentLocation loc : manager.getContentLocationsForDataBeanSaving(bean)) {
					msg += loc.getUrl() + " " + manager.getContentLength(loc) + " bytes, ";
				}
				throw new ContentLengthException(msg);
			} catch (IOException e1) {
				logger.error("another exception while handling " + Exceptions.getStackTrace(e), e);
			}					
		}
		bean.addContentLocation(location);		
	}
	

	public void removeContentLocationsFromDataBean(DataBean bean, StorageMethod method) {
		bean.removeContentLocations(method);
	}

	/**
	 * Get ContentLocations for DataBean. Only needed when saving a session.
	 * @param bean
	 */
	public List<ContentLocation> getContentLocationsForDataBeanSaving(DataBean bean) {
		return bean.getContentLocations();
	}

	
	public boolean uploadToStorageIfNeeded(DataBean bean) throws Exception {

		// check if already in storage
		if (Session.getSession().getServiceAccessor().getFileBrokerClient().isAvailable(bean.getId(), bean.getSize(), bean.getChecksum(), FileBrokerArea.STORAGE)) {
			return true;
		}
		
		// move from cache if possible
		if (Session.getSession().getServiceAccessor().getFileBrokerClient().moveFromCacheToStorage(bean.getId())) {
			return true;
		}
				
		// upload
		return upload(bean, FileBrokerArea.STORAGE, null);		
	}
	
	/**
	 * 
	 * @param bean
	 * @param progressListener
	 * @return 
	 * @throws NotEnoughDiskSpaceException
	 * @throws FileBrokerException
	 * @throws JMSException
	 * @throws IOException
	 * @throws Exception
	 */
	public boolean uploadToCacheIfNeeded(DataBean bean, CopyProgressListener progressListener) throws NotEnoughDiskSpaceException, FileBrokerException, JMSException, IOException, Exception {
		
		try {
			bean.getLock().readLock().lock();

			// upload only if not already available in cache or storage
			if (Session.getSession().getServiceAccessor().getFileBrokerClient().isAvailable(bean.getId(), bean.getSize(), bean.getChecksum(), FileBrokerArea.CACHE) ||
				Session.getSession().getServiceAccessor().getFileBrokerClient().isAvailable(bean.getId(), bean.getSize(), bean.getChecksum(), FileBrokerArea.STORAGE)) {
				return false;
			}
			
			// need to upload
			return upload(bean, FileBrokerArea.CACHE, progressListener);

		} finally {
			bean.getLock().readLock().unlock();
		}
	}

	private boolean upload(DataBean dataBean, FileBrokerArea area, CopyProgressListener progressListener) throws Exception {
		// check if content is still available
		if (dataBean.getContentLocations().size() == 0) {
			return false;
		}
		
		// try to upload
		try {
			String checksum = Session.getSession().getServiceAccessor().getFileBrokerClient().addFile(
					dataBean.getId(), 
					area, 
					getContentStream(dataBean, DataNotAvailableHandling.EXCEPTION_ON_NA), 
					getContentLength(dataBean), 
					progressListener);
			
			setOrVerifyChecksum(dataBean, checksum);

		} catch (Exception e) {
			logger.warn("could not upload data: " + dataBean.getName(), e);
			throw e;
		}
		return true;
	}

	/**
	 * Do the appropriate operation with checksum regardless of dataBean's state. When dataBean
	 * doesn't yet have checksum we can only set it. If checksum exists already, it must equal with 
	 * the given checksum or CheckusmException is thrown.
	 * 
	 * @param dataBean
	 * @param checksum
	 * @throws ChecksumException
	 */
	public void setOrVerifyChecksum(DataBean dataBean, String checksum)
			throws ChecksumException {
		if (checksum != null) {
			if (dataBean.getChecksum() == null) {
				dataBean.setChecksum(checksum);
			} else {
				if (!checksum.equals(dataBean.getChecksum())) {
					throw new ChecksumException();
				}
			}
		}
	}
	
	/**
	 * Do the appropriate operation with content length regardless of dataBean's state. When dataBean
	 * doesn't yet have content length we can only set it. If content length exists already, it must equal with 
	 * the given content length or ContentLengthException is thrown.
	 * 
	 * @param dataBean
	 * @param checksum
	 * @throws ContentLengthException
	 */
	public void setOrVerifyContentLength(DataBean dataBean, Long contentLength)
			throws ContentLengthException {
		if (contentLength != null) {
			if (dataBean.getSize() == null) {
				dataBean.setSize(contentLength);
			} else {
				if (!contentLength.equals(dataBean.getSize())) {
					throw new ContentLengthException();
				}
			}
		}
	}

	/**
	 * Returns the ContentLocation that is likely to be the fastest available and 
	 * supports random access. 
	 * All returned ContentLocations are checked to be accessible. Returns null
	 * if none of the locations are accessible.
	 */
	public ContentLocation getClosestRandomAccessContentLocation(DataBean bean) {

		List<ContentLocation> closestContentLocations = getClosestContentLocationList(bean);
			
		for (ContentLocation contentLocation : closestContentLocations) {
			if (contentLocation != null && contentLocation.method.isRandomAccess() && isAccessible(contentLocation)) {
				return contentLocation;
			}
		}

		// nothing was accessible
		return null;
	}
	
	/**
	 * Returns the ContentLocation that is likely to be the fastest available. 
	 * All returned ContentLocations are checked to be accessible. Returns null
	 * if none of the locations are accessible.
	 */
	private ContentLocation getClosestContentLocation(DataBean bean) {

		List<ContentLocation> closestContentLocations = getClosestContentLocationList(bean);
			
		for (ContentLocation contentLocation : closestContentLocations) {
			if (contentLocation != null && isAccessible(contentLocation)) {
				return contentLocation;
			}
		}

		// nothing was accessible
		return null;
	}

	
	private List<ContentLocation> getClosestContentLocationList(DataBean bean) {
		
		List<ContentLocation> closestLocations = new LinkedList<ContentLocation>();

		closestLocations.addAll(bean.getContentLocations(StorageMethod.LOCAL_FILE_METHODS));
		closestLocations.addAll(bean.getContentLocations(StorageMethod.REMOTE_FILE_METHODS));
		closestLocations.addAll(bean.getContentLocations(StorageMethod.OTHER_SLOW_METHODS));

		return closestLocations;
	}


	private boolean isAccessible(ContentLocation location) {
		return location.getHandler().isAccessible(location);
	}

	/**
	 * Get a approximate row count of files under MAX_BYTES_TO_COUNT in size.
	 * If there are more than MAX_ROWS_TO_COUNT, this number is returned. 
	 * 
	 * @param data
	 * @return
	 * @throws MicroarrayException
	 */
	public Long getFastRowCount(DataBean data) throws MicroarrayException {
		if (getContentLength(data) < MAX_BYTES_TO_COUNT && 			
				(data.hasTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES) || 
						data.hasTypeTag(BasicModule.TypeTags.TABLE_WITHOUT_COLUMN_NAMES) || 
						data.hasTypeTag(MicroarrayModule.TypeTags.PHENODATA))) {
			
			// check if rows are counted already 
			Long cachedCount = (Long)data.getFromContentBoundCache(AT_LEAST_ROWS_CACHENAME);
			if (cachedCount != null) {
				return cachedCount; 

			} else {
				// count rows
				try (Table rowCounter = data.queryFeatures("/column/*").asTable()) {
					long rowCount = 0;
					while (rowCounter != null && rowCounter.nextRow() && rowCount < MAX_ROWS_TO_COUNT) {
						rowCount++;
					}
					data.putToContentBoundCache(AT_LEAST_ROWS_CACHENAME, (Long)rowCount);
					return rowCount;
				}
			}			
		}
		return null;
	}
}
