package fi.csc.microarray.filebroker;

import java.io.File;

import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;


public class FileBrokerAreas {

	private File cacheRoot;
	private File storageRoot;

	public FileBrokerAreas(File repositoryRoot, String cachePath, String storagePath) {
		this.cacheRoot = new File(repositoryRoot, cachePath);
		this.storageRoot = new File(repositoryRoot, storagePath);
	}

	public boolean fileExists(String fileId, FileBrokerArea area) {

		// check id
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId)) {
			throw new IllegalArgumentException("illegal file id: " + fileId);
		}
		
		switch(area) {

		case CACHE:
			File cacheFile = new File(cacheRoot, fileId);
			if (cacheFile.exists()) {
				return true; 
			} else {
				return false;
			}
		case STORAGE:
			File storageFile = new File(storageRoot, fileId);
			if (storageFile.exists()) {
				return true; 
			} else {
				return false;
			}
		default:
			throw new IllegalArgumentException("illegal filebroker area: " + area);
		}
	}
	
	public boolean moveFromCacheToStorage(String fileId) {
		// check id
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId)) {
			throw new IllegalArgumentException("illegal file id: " + fileId);
		}

		File cacheFile = new File(cacheRoot, fileId);
		File storageFile = new File(storageRoot, fileId);
		
		return cacheFile.renameTo(storageFile);
	}

}
