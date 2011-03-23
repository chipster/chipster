package fi.csc.microarray.client.session;

import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
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

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.session.schema.FolderType;
import fi.csc.microarray.client.session.schema.ObjectFactory;
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

	private final int DATA_BLOCK_SIZE = 2048;
	
	private File sessionFile;
	private HashMap<DataBean, URL> newURLs;
	private int entryCounter = 0;

	private int itemIdCounter = 0;
	private HashMap<Integer, DataItem> itemIdMap = new HashMap<Integer, DataItem>();
	private HashMap<DataItem, Integer> reversedItemIdMap = new HashMap<DataItem, Integer>();

	private int operationIdCounter = 0;
	private HashMap<Integer, Operation> operationIdMap = new HashMap<Integer, Operation>();
	private HashMap<Operation, Integer> reversedOperationIdMap = new HashMap<Operation, Integer>();
	
	private DataManager dataManager = Session.getSession().getDataManager();

	ObjectFactory factory;
	SessionType sessionType;

	
	public SessionSaver(File sessionFile) {
		this.sessionFile = sessionFile;
		this.dataManager = Session.getSession().getDataManager();

	}
	
	public void saveSnapshot() throws IOException, JAXBException {

		this.newURLs = new HashMap<DataBean, URL>();
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

			// write data and gather metadata simultanously
			sessionType.setFormatVersion(ClientSession.SESSION_VERSION);

			// generate all ids
			generateIdsRecursively(dataManager.getRootFolder());

			// 1st pass, write most metadata
			saveRecursively(dataManager.getRootFolder(), zipOutputStream);

			// 2nd pass for links (if written in one pass, input dependent operation parameters break when reading)
//			saveLinksRecursively(dataManager.getRootFolder(), metadata);

			
			// validate meta data 
			
			// save meta data
			JAXBContext jaxbContext = JAXBContext.newInstance("fi.csc.microarray.client.session.schema");
			Marshaller marshaller = jaxbContext.createMarshaller();
			marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, true);
			marshaller.marshal(factory.createSession(sessionType), System.out);

			ZipEntry cpZipEntry = new ZipEntry(ClientSession.SESSION_METADATA_FILENAME);
			zipOutputStream.putNextEntry(cpZipEntry);
			marshaller.marshal(factory.createSession(sessionType), zipOutputStream);

			zipOutputStream.closeEntry() ;							

//			writeFile(zipOutputStream, METADATA_FILENAME, 
//			new ByteArrayInputStream(metadata.toString().getBytes()));

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
		Integer id = itemIdCounter++;
		itemIdMap.put(id, data);
		reversedItemIdMap.put(data, id);
	}
	
	
	private void saveRecursively(DataFolder folder, ZipOutputStream cpZipOutputStream) throws IOException {
		
		String folderId = fetchId(folder);
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
			folderType.setParent(fetchId(folder.getParent()));
		}
		sessionType.getFolder().add(folderType);
		
		//		metadata.append("DATAFOLDER " + folderId + "\n");
//		saveDataItemMetadata(folder, folderId, metadata);
	}	
	
	private void saveDataBeanMetadata(DataBean bean, URL newURL, String folderId) {
//		String beanId = fetchId(bean);
//		
//		// for now all data content goes to session --> type is local session
//		metadata.append("DATABEAN " + beanId + " " + newURL + " " + StorageMethod.LOCAL_SESSION + " " + bean.getRepositoryName() + "\n");
//		
//		if (bean.getOperation() != null) {
//			Operation operation = bean.getOperation();
//			String operId;
//			
//			// write operation or lookup already written
//			if (!operationIdMap.containsValue(operation) ) {
//				
//				operId = generateId(operation);
//				metadata.append("OPERATION " + operId + " " + operation.getID() + "\n");
//
//				// some parameters need inputs at loading time => write these first
//				if (operation.getBindings() != null) {
//					String beanIds = "";
//					for (DataBinding binding : operation.getBindings()) {
//						beanIds += fetchId(binding.getData()) + " ";
//					}
//					metadata.append("INPUTS " + operId + " " + beanIds + "\n");
//				}
//
//				for (Parameter parameter : operation.getParameters()) {
//
//					// Write parameter only when value is not empty
//					if (parameter.getValue() != null && !parameter.getValue().equals("")) {	
//						metadata.append("OPERATION_PARAMETER " + operId + " " +
//								parameter.getID() + " " + parameter.getValueAsString() + "\n");
//					}
//				}
//
//
//			} else {
//				operId = reversedOperationIdMap.get(operation).toString();
//			}
//			
//			metadata.append("OUTPUT " + operId + " " + beanId + "\n");
//			
//		}
//		
//		if (bean.getNotes() != null) {
//			// remove newlines from notes, they break loading
//			metadata.append("NOTES " + beanId + " " + bean.getNotes().replace('\n', ' ') + "\n");
//		}
//		
//		if (bean.getCacheUrl() != null) {
//			metadata.append("CACHED_URL " + beanId + " " + bean.getCacheUrl() + "\n");			
//		}
//		
//		saveDataItemMetadata(bean, beanId, metadata);
	}

	private void saveDataItemMetadata(DataItem data, String folderId) {
//		metadata.append("NAME " + folderId + " " + data.getName() + "\n");
//		if (data.getParent() != null) {
//			metadata.append("CHILD " + folderId + " " + fetchId(data.getParent()) + "\n");
//		}
	}

	
	private String fetchId(DataItem item) {
		return reversedItemIdMap.get(item).toString();
	}

	
	private void saveLinksRecursively(DataFolder folder) {
//		
//		for (DataItem data : folder.getChildren()) {
//			
//			if (data instanceof DataFolder) {
//				saveLinksRecursively((DataFolder)data, metadata);
//				
//			} else {
//				DataBean bean = (DataBean)data; 
//				for (Link type : Link.values()) {
//					for (DataBean target : bean.getLinkTargets(type)) {
//						String beanId = fetchId(bean);
//						String targetId = fetchId(target);				
//						metadata.append("LINK " + type.name() + " " + beanId + " " + targetId + "\n");
//					}
//				}		
//			}
//		}		
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
		Integer id = operationIdCounter++;
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
