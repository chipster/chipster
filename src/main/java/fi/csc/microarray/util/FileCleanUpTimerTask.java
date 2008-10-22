package fi.csc.microarray.util;

import java.io.File;
import java.util.TimerTask;

import org.apache.log4j.Logger;

public class FileCleanUpTimerTask extends TimerTask {

	private File baseDir;
	private long cutoff;
	
	private static final Logger logger = Logger.getLogger(FileCleanUpTimerTask.class);
	
	/**
	 * Walks the baseDir recursively and deleted files and directories older than cutoff.
	 * 
	 * If a directory is old but contains files (which are not too old), it is not deleted.
	 * 
	 * @param baseDir
	 * @param cutoff milliseconds 
	 */
	public FileCleanUpTimerTask(File baseDir, long cutoff) {
		this.baseDir = baseDir;
		this.cutoff = cutoff;
	}
	
	@Override
	public void run() {
		
		try {
			Files.cleanOldFiles(baseDir, cutoff);
		} catch (Throwable t) {
			logger.warn("Got throwable when cleaning up old files.", t);
		}
	}

	
}
	
