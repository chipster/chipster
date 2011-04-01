package fi.csc.microarray.client.session;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;

import org.apache.log4j.Logger;
import org.xml.sax.SAXException;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.session.schema.DataType;
import fi.csc.microarray.client.session.schema.FolderType;
import fi.csc.microarray.client.session.schema.LinkType;
import fi.csc.microarray.client.session.schema.NameType;
import fi.csc.microarray.client.session.schema.ObjectFactory;
import fi.csc.microarray.client.session.schema.OperationType;
import fi.csc.microarray.client.session.schema.ParameterType;
import fi.csc.microarray.client.session.schema.SessionType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.StorageMethod;
import fi.csc.microarray.databeans.handlers.ZipDataBeanHandler;
import fi.csc.microarray.util.IOUtils;

public class SessionSaver {

	private static final Logger logger = Logger.getLogger(SessionSaver.class);

	
	private final int DATA_BLOCK_SIZE = 2048;
	
	private File sessionFile;
	private HashMap<DataBean, URL> newURLs = new HashMap<DataBean, URL>();

	private int entryCounter = 0;

	private int itemIdCounter = 0;
	private HashMap<String, DataItem> itemIdMap = new HashMap<String, DataItem>();
	private HashMap<DataItem, String> reversedItemIdMap = new HashMap<DataItem, String>();

	private int operationIdCounter = 0;
	private HashMap<String, Operation> operationIdMap = new HashMap<String, Operation>();
	private HashMap<Operation, String> reversedOperationIdMap = new HashMap<Operation, String>();
	private HashMap<String, OperationType> operationTypeMap = new HashMap<String, OperationType>();
	
	
	private DataManager dataManager = Session.getSession().getDataManager();

	ObjectFactory factory;
	SessionType sessionType;

	
	public SessionSaver(File sessionFile) {
		this.sessionFile = sessionFile;
		this.dataManager = Session.getSession().getDataManager();

	}
	
	public void saveSession() throws IOException, JAXBException, SAXException {

		// figure out the target file
		boolean replaceOldSession = sessionFile.exists();
		File newSessionFile;
		File backupFile = null;
		if (replaceOldSession) {
			// TODO maybe avoid overwriting existing temp file
			newSessionFile = new File(sessionFile.getAbsolutePath() + "-temp.cs");
			backupFile = new File(sessionFile.getAbsolutePath() + "-backup.cs");
		} else {
			newSessionFile = sessionFile;
		}

		
		this.factory = new ObjectFactory();
		this.sessionType = factory.createSessionType();
		
		ZipOutputStream zipOutputStream = null;
		boolean createdSuccessfully = false;
		try {
			
			zipOutputStream = new ZipOutputStream(new BufferedOutputStream(new FileOutputStream(newSessionFile)));
			zipOutputStream.setLevel(1); // quite slow with bigger values														

			// save session version
			sessionType.setFormatVersion(ClientSession.SESSION_VERSION);

			// generate all ids
			generateIdsRecursively(dataManager.getRootFolder());

			// gather session data and save actual data to the zip file
			saveRecursively(dataManager.getRootFolder(), zipOutputStream);
			
			// validate and save session data
			Marshaller marshaller = ClientSession.getJAXBContext().createMarshaller();
			marshaller.setSchema(ClientSession.getSchema());
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
			marshaller.marshal(factory.createSession(sessionType), System.out);

			ZipEntry sessionDataZipEntry = new ZipEntry(ClientSession.SESSION_DATA_FILENAME);
			zipOutputStream.putNextEntry(sessionDataZipEntry);
			marshaller.marshal(factory.createSession(sessionType), zipOutputStream);

			zipOutputStream.closeEntry() ;							

			// FIXME finally?
			zipOutputStream.finish();
			zipOutputStream.close();

			
			
			// rename new session if replacing existing
			if (replaceOldSession) {
				
				// original to backup
				if (!sessionFile.renameTo(backupFile)) {
					throw new IOException("Creating backup " + sessionFile + " -> " + backupFile + " failed.");
				}
					
				// new to original
				if (newSessionFile.renameTo(sessionFile)) {
					createdSuccessfully = true;

					// remove backup
					backupFile.delete();
				} else {
					// try to move backup back to original
					// TODO remove new session file?
					if (backupFile.renameTo(sessionFile)) {
						throw new IOException("Moving new " + newSessionFile + " -> " + sessionFile + " failed, " +
								"restored original session file.");
					} else {
						throw new IOException("Moving new " + newSessionFile + " -> " + sessionFile + " failed, " +
						"also restoring original file failed, backup of original is " + backupFile);
					}
				}
			} 
			
			// session file is now saved, update the urls and handlers in the client
			for (DataBean bean: newURLs.keySet()) {

				// set new url and handler and type
				bean.setStorageMethod(StorageMethod.LOCAL_SESSION);
				bean.setContentUrl(newURLs.get(bean));
				bean.setHandler(new ZipDataBeanHandler());
			}

			createdSuccessfully = true;
			
		} catch (RuntimeException e) {
			// createdSuccesfully is false, so file will be deleted in finally block
			throw e;
			
		} catch (IOException e) {
			// createdSuccesfully is false, so file will be deleted in finally block
			throw e;

		} catch (JAXBException e) {
			// createdSuccesfully is false, so file will be deleted in finally block
			throw e;
		} finally {
			IOUtils.closeIfPossible(zipOutputStream); // called twice for normal execution, not a problem
			if (!replaceOldSession && !createdSuccessfully) {
				newSessionFile.delete(); // do not leave bad session files hanging around
			}
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
	
	
	private void saveRecursively(DataFolder folder, ZipOutputStream cpZipOutputStream) throws IOException {
		
		String folderId = reversedItemIdMap.get(folder);
		saveDataFolderMetadata(folder, folderId);
		
		for (DataItem data : folder.getChildren()) {
			if (data instanceof DataFolder) {
				saveRecursively((DataFolder)data, cpZipOutputStream);
				
			} else {
				DataBean bean = (DataBean)data;

				// create the new URL TODO check the ref
				String entryName = getNewEntryName();
				URL newURL = new URL(sessionFile.toURI().toURL(), "#" + entryName);

				// store the new URL temporarily
				newURLs.put(bean, newURL);

				// store metadata
				saveDataBeanMetadata(bean, newURL, folderId);
				
				// write bean contents to zip
				writeFile(cpZipOutputStream, entryName, bean.getContentByteStream());

			}
		}
	}


	private void saveDataFolderMetadata(DataFolder folder, String folderId) {
		FolderType folderType = factory.createFolderType();
		folderType.setId(folderId);
		folderType.setName(folder.getName());
		if (folder.getParent() != null) {
			String parentId = reversedItemIdMap.get(folder.getParent());
			if (parentId != null) {
				folderType.setParent(parentId);
			} else {
				logger.warn("unknown parent");
			}
		}
		
		if (folder.getChildCount() > 0) {
			for (DataItem child : folder.getChildren()) {
				String childId = reversedItemIdMap.get(child);
				if (childId != null) { 
					folderType.getChild().add(reversedItemIdMap.get(child));
				} else {
					logger.warn("unknown child: " + child.getName());
				}
			}
		}
		
		sessionType.getFolder().add(folderType);
	}	
	
	private void saveDataBeanMetadata(DataBean bean, URL newURL, String folderId) {
		String beanId = reversedItemIdMap.get(bean);
		
		// save the basic data
		DataType dataType = factory.createDataType();
		dataType.setId(beanId);
		dataType.setName(bean.getName());
		if (bean.getParent() != null) {
			String parentId = reversedItemIdMap.get(bean.getParent());
			if (parentId != null) {
				dataType.setFolder(parentId);
			} else {
				logger.warn("unknown parent");
			}
		}
		
		dataType.setNotes(bean.getNotes());
		
		// storage method
		// for now all data content goes to session --> type is local session
		dataType.setStorageType(StorageMethod.LOCAL_SESSION.name());
		
		// url
		dataType.setUrl(newURL.toString());
		
		// cache url
		if (bean.getCacheUrl() != null) {
			dataType.setCacheUrl(bean.getCacheUrl().toString());
		}
		

		// FIXME accept beans without operation?
		if (bean.getOperation() != null) {
			Operation operation = bean.getOperation();
			String operId;
			
			// write operation or lookup already written
			if (!operationIdMap.containsValue(operation) ) {
				operId = generateId(operation);
				saveOperationMetadata(operation, operId);

			} else {
				operId = reversedOperationIdMap.get(operation).toString();
			}

			// link data to operation
			operationTypeMap.get(operId).getOutput().add(beanId);
			
			// link the operation to data
			dataType.setResultOf(operId);
			
		}
		
		// links to other datasets
		for (Link type : Link.values()) {
			for (DataBean target : bean.getLinkTargets(type)) {
				// FIXME check for targetId not null
				String targetId = reversedItemIdMap.get(target);				
				LinkType linkType = factory.createLinkType();
				linkType.setTarget(targetId);
				linkType.setType(type.name());
				
				dataType.getLink().add(linkType);
			}
		}		
		

		sessionType.getData().add(dataType);

	}

	
	private void saveOperationMetadata(Operation operation, String operationId) {
		OperationType operationType = factory.createOperationType();
		
		// session id
		operationType.setId(operationId);
		
		// name
		NameType nameType = factory.createNameType();
		nameType.setId(operation.getID());
		nameType.setDisplayName(operation.getDisplayName());
		operationType.setName(nameType);
		
		// parameters
		for (Parameter parameter : operation.getParameters()) {

			// Write parameter only when value is not empty
			if (parameter.getValue() != null && !parameter.getValue().equals("")) {	
				ParameterType parameterType = factory.createParameterType();
				NameType parameterNameType = factory.createNameType();
				parameterNameType.setId(parameter.getID());
				parameterNameType.setDisplayName(parameter.getDisplayName());
				parameterType.setName(parameterNameType);
				parameterType.setValue(parameter.getValueAsString());
				operationType.getParameter().add(parameterType);
			}
		}

		// inputs
		for (DataBinding binding : operation.getBindings()) {
			// FIXME check inputId for null
			String inputId = reversedItemIdMap.get(binding.getData());
			operationType.getInput().add(inputId); 
		}
		
		sessionType.getOperation().add(operationType);
		operationTypeMap.put(operationId, operationType);
	}	


	private String getNewEntryName() {
		return "file-" + entryCounter++;
	}

	private void writeFile(ZipOutputStream out, String name, InputStream in) throws IOException {
		
		int byteCount;
		ZipEntry cpZipEntry = new ZipEntry(name);
		out.putNextEntry(cpZipEntry);

		byte[] b = new byte[DATA_BLOCK_SIZE];

		while ( (byteCount = in.read(b, 0, DATA_BLOCK_SIZE)) != -1 ) {
			out.write(b, 0, byteCount);
		}

		out.closeEntry() ;							
	}

	private String generateId(Operation operation) {
		String id = String.valueOf(operationIdCounter);
		operationIdCounter++;
		operationIdMap.put(id, operation);
		reversedOperationIdMap.put(operation, id);
		return id.toString();
	}

	
	public static void main(String[] args) throws JAXBException {
		ObjectFactory objFactory = new ObjectFactory();
		SessionType session = objFactory.createSessionType();
		session.setFormatVersion(3);
		FolderType folder = objFactory.createFolderType();
		folder.setName("Jou folder");
		folder.setId("folder-1");
		session.getFolder().add(folder);
		
		JAXBContext jaxbContext = JAXBContext.newInstance("fi.csc.microarray.client.session.schema");
		Marshaller marshaller = jaxbContext.createMarshaller();
		marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
		marshaller.marshal(objFactory.createSession(session), System.out);

	}
}
