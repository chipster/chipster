package fi.csc.microarray.client.session;

import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map.Entry;

import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import org.apache.log4j.Logger;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import de.schlichtherle.truezip.zip.ZipEntry;
import de.schlichtherle.truezip.zip.ZipOutputStream;
import fi.csc.microarray.client.NameID;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.InputRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.session.schema2.DataType;
import fi.csc.microarray.client.session.schema2.FolderType;
import fi.csc.microarray.client.session.schema2.InputType;
import fi.csc.microarray.client.session.schema2.LinkType;
import fi.csc.microarray.client.session.schema2.LocationType;
import fi.csc.microarray.client.session.schema2.NameType;
import fi.csc.microarray.client.session.schema2.ObjectFactory;
import fi.csc.microarray.client.session.schema2.OperationType;
import fi.csc.microarray.client.session.schema2.ParameterType;
import fi.csc.microarray.client.session.schema2.SessionType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.ContentLocation;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataManager.StorageMethod;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.SwingTools;
import fi.csc.microarray.util.UrlTransferUtil;

/**
 * @author hupponen
 *
 */
public class SessionSaver {

	private static final Logger logger = Logger.getLogger(SessionSaver.class);

	
	private final int DATA_BLOCK_SIZE = 2048;
	
	private File sessionFile;
	private URL sessionUrl;
	private HashMap<DataBean, URL> newURLs = new HashMap<DataBean, URL>();

	private int entryCounter = 0;
	private int sourceCodeEntryCounter = 0;
	
	private int itemIdCounter = 0;
	private HashMap<String, DataItem> itemIdMap = new HashMap<String, DataItem>();
	private HashMap<DataItem, String> reversedItemIdMap = new HashMap<DataItem, String>();

	private int operationIdCounter = 0;
	private HashMap<String, OperationRecord> operationRecordIdMap = new HashMap<String, OperationRecord>();
	private HashMap<OperationRecord, String> reversedOperationRecordIdMap = new HashMap<OperationRecord, String>();
	private HashMap<String, OperationType> operationRecordTypeMap = new HashMap<String, OperationType>();
	
	private DataManager dataManager;

	private ObjectFactory factory;
	private SessionType sessionType;

	private String validationErrors;



	/**
	 * Create a new instance for every session to be saved.
	 * 
	 * @param sessionFile file to write out metadata and possible data
	 */
	public SessionSaver(File sessionFile, DataManager dataManager) {
		this.sessionFile = sessionFile;
		this.sessionUrl = null;
		this.dataManager = dataManager;

	}

	/**
	 * Create a new instance for every session to be saved.
	 * 
	 * @param sessionUrl url to write out metadata
	 */
	public SessionSaver(URL sessionUrl, DataManager dataManager) {
		this.sessionFile = null;
		this.sessionUrl = sessionUrl;
		this.dataManager = dataManager;

	}

	/**
	 * Use getValidationException() to get the reason for failed validation
	 * of the metadata (when returning false).
	 * 
	 * @return true if the written metadata was valid
	 * @throws Exception if something else than validation fails
	 */
	public boolean saveSession() throws Exception{

		gatherMetadata(true, true);
		boolean metadataValid = validateMetadata();
	
		writeSessionToFile(true);
		updateDataBeanURLsAndHandlers();
		
		return metadataValid;
	}

	public void saveLightweightSession() throws Exception {

		gatherMetadata(false, false);
		writeSessionToFile(false);
	}

	public LinkedList<URL> saveStorageSession() throws Exception {

		// move data bean contents to storage
		LinkedList<URL> urls = new LinkedList<URL>();
		for (DataBean dataBean : dataManager.databeans()) {
			URL url = dataManager.putToStorage(dataBean);
			if (url != null) {
				urls.add(url);				
			}
		}
		
		// save metadata
		gatherMetadata(false, true);
		writeSessionToUrl(false);
		
		return urls;
	}
	
	
	/**
	 * Gather the metadata form the data beans, folders and operations.
	 * 
	 * @throws IOException
	 * @throws JAXBException
	 */
	private void gatherMetadata(boolean saveData, boolean skipLocalLocations) throws IOException, JAXBException {
		// xml schema object factory and xml root
		this.factory = new ObjectFactory();
		this.sessionType = factory.createSessionType();

		// save session version
		sessionType.setFormatVersion(UserSession.SESSION_VERSION);

		// generate all ids
		generateIdsRecursively(dataManager.getRootFolder());

		// gather meta data
		saveMetadataRecursively(dataManager.getRootFolder(), saveData, skipLocalLocations);
	}


	/**
	 * 
	 * @throws JAXBException
	 * @throws SAXException
	 */
	private boolean validateMetadata() throws JAXBException, SAXException {
		Marshaller marshaller = UserSession.getJAXBContext().createMarshaller();
		marshaller.setSchema(UserSession.getSchema());
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
		
		NonStoppingValidationEventHandler validationEventHandler = new NonStoppingValidationEventHandler();
		marshaller.setEventHandler(validationEventHandler);
		
		marshaller.marshal(factory.createSession(sessionType), new DefaultHandler());
	
		if (!validationEventHandler.hasEvents()) {
			 return true;
		} else {
			this.validationErrors = validationEventHandler.getValidationEventsAsString();
			return false;
		}
		 
	}

	/**
	 * Write metadata over URL. 
	 * 
	 * @param saveData if true, also actual contents of databeans are saved 
	 * 
	 */
	private void writeSessionToUrl(boolean saveData) throws Exception {

		HttpURLConnection conn = UrlTransferUtil.prepareForUpload(sessionUrl);
		
		try {
			OutputStream out = conn.getOutputStream();
			try {
				writeSessionContents(saveData, out); 
			} finally {
				IOUtils.closeIfPossible(out);
			}
			
			// need to check this to guarantee upload
			if (!UrlTransferUtil.isSuccessfulCode(conn.getResponseCode())) {
    			throw new IOException("PUT was not successful: " + conn.getResponseCode() + " " + conn.getResponseMessage());
    		}
			
		} finally {
			IOUtils.disconnectIfPossible(conn); 
		}
		
	}

	/**
	 * Write the metadata file and possibly data bean contents to the zip file.
	 * 
	 * @param saveData if true, also actual contents of databeans are saved 
	 * 
	 */
	private void writeSessionToFile(boolean saveData) throws Exception {

		// figure out the target file, use temporary file if target already exists
		boolean replaceOldSession = sessionFile.exists();
		File newSessionFile;
		File backupFile = null;
		if (replaceOldSession) {
			newSessionFile = new File(sessionFile.getAbsolutePath() + "-save-temp.zip");
			backupFile = new File(sessionFile.getAbsolutePath() + "-save-backup.zip");
		} else {
			newSessionFile = sessionFile;
		}

		// write data to zip file 
		FileOutputStream out = null;
		try {
			out = new FileOutputStream(newSessionFile);
			writeSessionContents(saveData, out);
			IOUtils.closeIfPossible(out);
			
		} catch (Exception e) {
			IOUtils.closeIfPossible(out);
			// don't leave the new session file lying around if something went wrong
			newSessionFile.delete();
			throw e;
		}
		
		
		// rename new session if replacing existing
		if (replaceOldSession) {

			// close open files, because they might stop us from overwriting existing session file
			dataManager.flushSession();
			
			// original to backup
			if (!sessionFile.renameTo(backupFile)) {
				throw new IOException("Creating backup file " + backupFile.getAbsolutePath() + " failed.");
			}

			// new to original
			if (newSessionFile.renameTo(sessionFile)) {

				// remove backup
				backupFile.delete();
			} else {
				// try to move backup back to original
				if (backupFile.renameTo(sessionFile)) {
					throw new IOException("Moving new session file " + newSessionFile + " -> " + sessionFile + " failed, " +
					"restored original session file.");
				} else {
					throw new IOException("Moving new session file " + newSessionFile + " -> " + sessionFile + " failed, " +
							"also restoring original file failed, backup of original is " + backupFile);
				}
			}
		} 
	}

	private void writeSessionContents(boolean saveData, OutputStream out) throws Exception {

		ZipOutputStream zipOutputStream = null;
		try {	
			zipOutputStream = new ZipOutputStream(new BufferedOutputStream(out));
			zipOutputStream.setLevel(1); // quite slow with bigger values														

			// save meta data
			ZipEntry sessionDataZipEntry = new ZipEntry(UserSession.SESSION_DATA_FILENAME);
			zipOutputStream.putNextEntry(sessionDataZipEntry);
			Marshaller marshaller = UserSession.getJAXBContext().createMarshaller();
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
			// TODO disable validation
			marshaller.setEventHandler(new NonStoppingValidationEventHandler());
			marshaller.marshal(factory.createSession(sessionType), zipOutputStream);
			zipOutputStream.closeEntry() ;							

			// save data bean contents
			if (saveData) {
				writeDataBeanContentsToZipFile(zipOutputStream);
			}
			
			// save source codes
			writeSourceCodesToZip(zipOutputStream);
			
			// close the zip stream
			zipOutputStream.close();
		} 
		
		catch (Exception e) {
			IOUtils.closeIfPossible(zipOutputStream);
			throw e;
		}
	}


	/**
	 * After the session file has been saved, update the urls and handlers in the client
	 * to point to the data inside the session file.
	 *
	 */
	private void updateDataBeanURLsAndHandlers() {
		for (DataBean bean: newURLs.keySet()) {
			// set new url and handler and type
			dataManager.addUrl(bean, StorageMethod.LOCAL_SESSION, newURLs.get(bean));
		}
	}
	
	
	private int generateIdsRecursively(DataFolder folder) throws IOException {
		
		int dataCount = 0;
		
		generateId(folder);
		
		for (DataItem data : folder.getChildren()) {
			if (data instanceof DataFolder) {
				int recDataCount = generateIdsRecursively((DataFolder)data);
				dataCount += recDataCount;
				
			} else {
				generateId((DataBean)data);
				dataCount++;
			}
		}

		return dataCount;
	}

	private void generateId(DataItem data) {
		String id = String.valueOf(itemIdCounter);
		itemIdCounter++;
		itemIdMap.put(id, data);
		reversedItemIdMap.put(data, id);
	}
	
	
	private String generateId(OperationRecord operationRecord) {
		String id = String.valueOf(operationIdCounter);
		operationIdCounter++;
		operationRecordIdMap.put(id, operationRecord);
		reversedOperationRecordIdMap.put(operationRecord, id);
		return id.toString();
	}

	private void saveMetadataRecursively(DataFolder folder, boolean saveData, boolean skipLocalLocations) throws IOException {
		
		String folderId = reversedItemIdMap.get(folder);
		saveDataFolderMetadata(folder, folderId);
		
		for (DataItem data : folder.getChildren()) {
			if (data instanceof DataFolder) {
				saveMetadataRecursively((DataFolder)data, saveData, skipLocalLocations);
				
			} else {
				DataBean bean = (DataBean)data;

				// create the new URL
				String entryName = getNewZipEntryName();
				URL newURL = null;
				
				if (saveData) {

					// data is saved to zip, change URL to point there 
					newURL = new URL(sessionFile.toURI().toURL(), "#" + entryName);

					// store the new URL temporarily
					newURLs.put(bean, newURL);
				}
				
				// store metadata
				saveDataBeanMetadata(bean, newURL, folderId, skipLocalLocations);

			}
		}
	}


	private void saveDataFolderMetadata(DataFolder folder, String folderId) {
		FolderType folderType = factory.createFolderType();
		
		// name
		folderType.setId(folderId);
		folderType.setName(folder.getName());
		
		// parent
		if (folder.getParent() != null) {
			String parentId = reversedItemIdMap.get(folder.getParent());
			if (parentId != null) {
				folderType.setParent(parentId);
			} else {
				logger.warn("unknown parent");
			}
		}
		
		// children
		if (folder.getChildCount() > 0) {
			for (DataItem child : folder.getChildren()) {
				String childId = reversedItemIdMap.get(child);
				if (childId != null) { 
					folderType.getChild().add(childId);
				} else {
					logger.warn("unknown child: " + child.getName());
				}
			}
		}
		
		sessionType.getFolder().add(folderType);
	}	
	
	
	private void saveDataBeanMetadata(DataBean bean, URL newURL, String folderId, boolean skipLocalLocations) {
		String beanId = reversedItemIdMap.get(bean);
		DataType dataType = factory.createDataType();
	
		// name and id
		dataType.setId(beanId);
		dataType.setName(bean.getName());

		// parent
		if (bean.getParent() != null) {
			String parentId = reversedItemIdMap.get(bean.getParent());
			if (parentId != null) {
				dataType.setFolder(parentId);
			} else {
				logger.warn("unknown parent");
			}
		}
		
		// notes
		dataType.setNotes(bean.getNotes());

		// write all URL's
		for (ContentLocation location : bean.getContentLocations()) {
			if (skipLocalLocations && location.getMethod().isLocal()) {
				continue; // do not save local locations to remote sessions
			}
			LocationType locationType = new LocationType();
			locationType.setMethod(location.getMethod().toString());
			locationType.setUrl(location.getUrl().toString());
			dataType.getLocation().add(locationType);
		}
		
		// write newly created URL inside session files, if exists
		if (newURL != null) {
			LocationType locationType = new LocationType();
			locationType.setMethod(StorageMethod.LOCAL_SESSION.name());
			locationType.setUrl("file:#" + newURL.getRef());
			dataType.getLocation().add(locationType);
		}
		
		// for now, accept beans without operation
		if (bean.getOperationRecord() != null) {
			OperationRecord operationRecord = bean.getOperationRecord();
			String operId;
			
			// write operation or lookup already written
			if (!operationRecordIdMap.containsValue(operationRecord) ) {
				operId = generateId(operationRecord);
				saveOperationMetadata(operationRecord, operId);

			} else {
				operId = reversedOperationRecordIdMap.get(operationRecord).toString();
			}

			// link data to operation
			operationRecordTypeMap.get(operId).getOutput().add(beanId);
			
			// link the operation to data
			dataType.setResultOf(operId);
		}
		
		// links to other datasets
		for (Link type : Link.values()) {
			for (DataBean target : bean.getLinkTargets(type)) {
				String targetId = reversedItemIdMap.get(target);				
				// if for some weird reason target was not around when generating ids, skip it
				if (targetId == null) {
					continue;
				}
				LinkType linkType = factory.createLinkType();
				linkType.setTarget(targetId);
				linkType.setType(type.name());
				
				dataType.getLink().add(linkType);
			}
		}		
		
		sessionType.getData().add(dataType);
	}

	
	private void saveOperationMetadata(OperationRecord operationRecord, String operationId) {
		OperationType operationType = factory.createOperationType();
		
		// session id
		operationType.setId(operationId);
		
		// name
		NameType nameType = createNameType(operationRecord.getNameID());
		operationType.setName(nameType);
		
		// parameters
		for (ParameterRecord parameterRecord : operationRecord.getParameters()) {

			// Write parameter only when value is not empty
			if (parameterRecord.getValue() != null && !parameterRecord.getValue().equals("")) {	
				ParameterType parameterType = factory.createParameterType();
				parameterType.setName(createNameType(parameterRecord.getNameID()));
				parameterType.setValue(parameterRecord.getValue());
				operationType.getParameter().add(parameterType);
			}
		}

		// inputs
		for (InputRecord inputRecord : operationRecord.getInputs()) {

			String inputID = reversedItemIdMap.get(inputRecord.getValue());
			// skip inputs which were not around when generating ids
			if (inputID == null) {
				continue;
			}
			InputType inputType = factory.createInputType();
			inputType.setName(createNameType(inputRecord.getNameID()));
			inputType.setData(inputID);
			
			operationType.getInput().add(inputType);
		}
		
		// category
		operationType.setCategory(operationRecord.getCategoryName());
		if (operationRecord.getCategoryColor() != null) {
			operationType.setCategoryColor(SwingTools.colorToHexString(operationRecord.getCategoryColor()));
		}

		// module
		operationType.setModule(operationRecord.getModule());

		// source code
		if (operationRecord.getSourceCode() != null) {
			String entryName = getNewSourceCodeEntryName(operationRecord.getNameID().getID());
			operationType.setSourceCodeFile(entryName);
		}
		
		sessionType.getOperation().add(operationType);
		operationRecordTypeMap.put(operationId, operationType);
	}	


	private String getNewZipEntryName() {
		return "file-" + entryCounter++;
	}

	private String getNewSourceCodeEntryName(String prefix) {
		return "source-code-" + sourceCodeEntryCounter++ + "-" + prefix + ".txt";
	}
	
	private void writeDataBeanContentsToZipFile(ZipOutputStream zipOutputStream) throws IOException {
		for (Entry<DataBean, URL> entry : this.newURLs.entrySet()) {
			String entryName = entry.getValue().getRef();

			// write bean contents to zip
			writeFile(zipOutputStream, entryName, entry.getKey().getContentStream(DataNotAvailableHandling.EXCEPTION_ON_NA));
		}
	}
	
	private void writeSourceCodesToZip(ZipOutputStream zipOutputStream) throws IOException {
		for (Entry<String, OperationType> entry : this.operationRecordTypeMap.entrySet()) {
			String sourceCodeFileName = entry.getValue().getSourceCodeFile();
			if (sourceCodeFileName != null && !sourceCodeFileName.isEmpty()) {
				writeFile(zipOutputStream, sourceCodeFileName, new ByteArrayInputStream(this.operationRecordIdMap.get(entry.getKey()).getSourceCode().getBytes()));
			}
		}
	}
	
	private void writeFile(ZipOutputStream out, String name, InputStream in) throws IOException {
		
		int byteCount;
		ZipEntry cpZipEntry = new ZipEntry(name);
		out.putNextEntry(cpZipEntry);

		byte[] b = new byte[DATA_BLOCK_SIZE];

		while ( (byteCount = in.read(b, 0, DATA_BLOCK_SIZE)) != -1 ) {
			out.write(b, 0, byteCount);
		}

		out.closeEntry();							
	}

	
	private NameType createNameType(String id, String displayName, String desription) {
		NameType nameType = factory.createNameType();
		nameType.setId(id);
		nameType.setDisplayName(displayName);
		nameType.setDescription(desription);
		return nameType;
	}
	

	private NameType createNameType(NameID nameID) {
		return createNameType(nameID.getID(), nameID.getDisplayName(), nameID.getDescription());
	}

	public String getValidationErrors() {
		return this.validationErrors;
	}
	
	public static void dumpSession(DataFolder folder, StringBuffer buffer) {

		for (DataItem data : folder.getChildren()) {
			
			if (data instanceof DataFolder) {
				dumpSession((DataFolder)data, buffer);
				
			} else {
				DataBean bean = (DataBean)data;
				buffer.append("\nBean: " + bean.getName() + "\n");
				
				for (ContentLocation locations : bean.getContentLocations()) {
					buffer.append("  " + locations.getMethod() + ": \t" + locations.getUrl() + "\n");
				}

			}
		}
		
	}
	
//    public static void main(String args[])
//    {                
//            try
//            {
//                    String zipFile = "/home/akallio/Desktop/test.zip";
//                   
//                    //create byte buffer
//                    byte[] buffer = new byte[1024];
//                    /*
//                     * To create a zip file, use
//                     *
//                     * ZipOutputStream(OutputStream out)
//                     * constructor of ZipOutputStream class.
//                     *  
//                     */
//                     
//                     //create object of FileOutputStream
//                     FileOutputStream fout = new FileOutputStream(zipFile);
//                     
//                     //create object of ZipOutputStream from FileOutputStream
//                     ZipOutputStream zout = new ZipOutputStream(fout);
//                     
//                     /*
//                      * To begin writing ZipEntry in the zip file, use
//                      *
//                      * void putNextEntry(ZipEntry entry)
//                      * method of ZipOutputStream class.
//                      *
//                      * This method begins writing a new Zip entry to
//                      * the zip file and positions the stream to the start
//                      * of the entry data.
//                      */
//                     
//                     zout.putNextEntry(new ZipEntry("jee"));
//                     
//                     /*
//                      * After creating entry in the zip file, actually
//                      * write the file.
//                      */
//                     int length = buffer.length;
//                     
//                     zout.write(buffer, 0, length);
//                     
//                     /*
//                      * After writing the file to ZipOutputStream, use
//                      *
//                      * void closeEntry() method of ZipOutputStream class to
//                      * close the current entry and position the stream to
//                      * write the next entry.
//                      */
//                     
//                      zout.closeEntry();
//                     
//                      //close the ZipOutputStream
//                      zout.close();
//                     
//                      System.out.println("Zip file has been created!");
//           
//            }
//            catch(IOException ioe)
//            {
//                    System.out.println("IOException :" + ioe);
//            }
//           
//    }
}
