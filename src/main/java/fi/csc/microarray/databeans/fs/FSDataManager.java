package fi.csc.microarray.databeans.fs;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;

import org.mortbay.util.IO;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataItemCreatedEvent;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataManagerBase;
import fi.csc.microarray.databeans.DataBean.Link;

public class FSDataManager extends DataManagerBase {

	private FSDataFolder rootFolder;	
	private File repositoryRoot;
	
	public FSDataManager() throws IOException {
		rootFolder = (FSDataFolder)createFolder(DataManager.ROOT_NAME);

		// initialize repository 		
		repositoryRoot = createRepository();
	}
	

	public void setRootFolder(FSDataFolder folder) {
		this.rootFolder = folder;		
	}
	
	public DataFolder getRootFolder() {
		return rootFolder;
	}

	public DataFolder createFolder(String name) {
		DataFolder folder = new FSDataFolder(this, name);
		return folder;
	}

	public DataFolder createFolder(DataFolder root, String name) {
		DataFolder folder = new FSDataFolder(this, name);
		root.addChild(folder); // events are dispatched from here
		return folder;
	}


	public DataBean createDataBean(String name) throws MicroarrayException {
		return createDataBean(name, null, new DataBean[] {});
	}

	public DataBean createDataBean(String name, DataFolder folder, DataBean... sources) throws MicroarrayException {
		File contentFile;
		try {
			contentFile = createNewRepositoryFile(name);
		} catch (IOException e) {
			throw new MicroarrayException(e);
		}
		
		return createDataBean(name, folder, sources, contentFile);
	}
	
	public FSDataBean createDataBean(String name, InputStream content) throws MicroarrayException {
		return createDataBean(name, content, null, new DataBean[] {});
	}
	
	public FSDataBean createDataBean(String name, InputStream content, DataFolder folder, DataBean... sources) throws MicroarrayException {

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
		FSDataBean bean = new FSDataBean(name, guessContentType(name), new Date(), sources, folder, this, contentFile);

		dispatchEventIfVisible(new DataItemCreatedEvent(bean));
		
		return bean;
	}

	/**
	 * Use this only if you have first created the new file with the
	 * createNewRepositoryFile(String name) method.
	 * 
	 */
	public FSDataBean createDataBean(String name, File contentFile) throws MicroarrayException {		
		return createDataBean(name, null, new DataBean[] {}, contentFile);
	}

	/**
	 * Use this only if you have first created the new file with the
	 * createNewRepositoryFile(String name) method.
	 * 
	 */
	public FSDataBean createDataBean(String name, DataFolder folder, DataBean[] sources, File contentFile) throws MicroarrayException {

		FSDataBean dataBean = new FSDataBean(name, guessContentType(name), new Date(), sources, folder, this, contentFile);
		dispatchEventIfVisible(new DataItemCreatedEvent(dataBean));
		return dataBean;
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
		// check the file name 
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


	public List<DataItem> loadSnapshot(File snapshotDir, DataFolder parentFolder, ClientApplication application) throws IOException, MicroarrayException {
		FSSnapshottingSession session = new FSSnapshottingSession(this, application);
		List<DataItem> newItems = session.loadFromSnapshot(snapshotDir, parentFolder);
		return newItems;
	}

	public List<DataItem> loadOldSnapshot(File snapshotDir, DataFolder parentFolder, ClientApplication application) throws IOException, MicroarrayException {
		OldFSSnapshottingSession session = new OldFSSnapshottingSession(this);
		List<DataItem> newItems = session.loadFromSnapshot(snapshotDir, parentFolder);
		return newItems;
	}


	public int saveSnapshot(File snapshotDir, ClientApplication application) throws IOException {
		FSSnapshottingSession session = new FSSnapshottingSession(this, application);
		return session.saveSnapshot(snapshotDir);
	}


	public void delete(DataItem data) {
		
		if (data instanceof DataFolder) {
			deleteDataFolder((DataFolder)data);
			
		} else {
			deleteDataBean((DataBean)data);
		}		
	}
	
	private void deleteDataBean(DataBean bean) {

		FSDataBean fsDataBean = (FSDataBean)bean;
		
		// remove links
		for (Link linkType : Link.values()) {
			// Remove outgoing links
			for (DataBean target : fsDataBean.getLinkTargets(linkType)) {
				fsDataBean.removeLink(linkType, target);
			}
			// Remove incoming links
			for (DataBean source : fsDataBean.getLinkSources(linkType)) {
				source.removeLink(linkType, fsDataBean);
			}
		}

		// remove this bean
		DataFolder folder = fsDataBean.getParent();
		if (folder != null) {
			folder.removeChild(fsDataBean);
		}
		
		// remove physical file
		fsDataBean.delete();
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

}