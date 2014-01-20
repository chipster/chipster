package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;

import org.apache.commons.io.FileUtils;

import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;


public class FileBrokerAreas {

	private File cacheRoot;
	private File storageRoot;

	public FileBrokerAreas(File repositoryRoot, String cachePath, String storagePath) {
		this.cacheRoot = new File(repositoryRoot, cachePath);
		this.storageRoot = new File(repositoryRoot, storagePath);
	}
	
	public Long getSize(String fileId, FileBrokerArea area) {
		if (fileExists(fileId, area)) {
			return getFile(fileId, area).length();
		}
		return null;
	}

	public boolean fileExists(String fileId, FileBrokerArea area) {
		return getFile(fileId, area).exists();
	}
	
	private File getFile(String fileId, FileBrokerArea area) {

		// check id
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId)) {
			throw new IllegalArgumentException("illegal file id: " + fileId);
		}
		
		switch(area) {

		case CACHE:
			File cacheFile = new File(cacheRoot, fileId);
			return cacheFile;
		case STORAGE:
			File storageFile = new File(storageRoot, fileId);
			return storageFile;
		default:
			throw new IllegalArgumentException("illegal filebroker area: " + area);
		}
	}
	
	public boolean moveFromCacheToStorage(String fileId) throws IOException {
		// check id
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId)) {
			throw new IllegalArgumentException("illegal file id: " + fileId);
		}

		File cacheFile = new File(cacheRoot, fileId);
		File storageFile = new File(storageRoot, fileId);
		
		return renameDataFile(cacheFile, storageFile);
	}
	
	/**
	 * Rename or move data file and its md5 file. Ensures that md5 file is accessible
	 * all the time during the operation.  
	 * 
	 * @param sourceDataFile
	 * @param destDataFile
	 * @return
	 * @throws IOException
	 */
	private boolean renameDataFile(File sourceDataFile, File destDataFile) throws IOException {
		
		
		File sourceMd5File = Md5FileUtils.getMd5File(sourceDataFile);
		File destMd5File = Md5FileUtils.getMd5File(destDataFile);
		
		boolean sameDirectory = sourceDataFile.getParentFile().getAbsoluteFile().equals(destDataFile.getParentFile().getAbsoluteFile());
		
		// don't try to copy if md5 checksums are disabled
		boolean copyChecksum = sourceMd5File.exists() && !sameDirectory;

		if (copyChecksum) {
			FileUtils.copyFile(sourceMd5File, destMd5File);
		}
		
		boolean renameSuccess = sourceDataFile.renameTo(destDataFile);
		
		if (copyChecksum) {
			if (renameSuccess) {
				sourceMd5File.delete();			
			} else {
				destMd5File.delete();
			}
		}
		
		return renameSuccess;
	}

	public String getChecksum(String fileId, FileBrokerArea area) throws ChecksumParseException, IOException {
		return Md5FileUtils.readMd5(getFile(fileId, area));
	}

}
