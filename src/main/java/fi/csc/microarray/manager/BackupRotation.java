package fi.csc.microarray.manager;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;

import org.apache.log4j.Logger;
import org.joda.time.DateTime;
import org.joda.time.format.ISODateTimeFormat;

/**
 * Deletion logic for backup files. There are multiple overlapping rules for files that are kept:
 * <ul>
 * <li>last 30 files</li>
 * <li>files of last 30 days</li>
 * <li>first file of each month</li>
 * </ul>
 * 
 * With daily backups the first two rules are pretty much the same, but offer some additional protection
 * against misconfiguration of the backup timer or server clock. First file of each month is stored forever and
 * thus not really rotated. 
 * 
 * @author klemela
 */
public class BackupRotation {
	
	private final Logger logger = Logger.getLogger(BackupRotation.class);

	private File backupRoot;
	
	private static final int ROTATION_COUNT = 30;
	private static final int ROTATION_DAYS = 30;
	
	private static final String filenamePrefix = "chipster-manager-db-backup-";
	private static final String filenamePostfix = ".zip";

	public BackupRotation(File backupRoot) {
		this.backupRoot = backupRoot;
	}

	/**
	 * Last 30 files, files newer than 30 days and first file of each month is kept, all other files
	 * are deleted (provided we can parse the filename).
	 * 
	 * @return number of files that were tried to delete 
	 */
	public int rotate() {
		/* 
		 * The logic here is a little bit deceptive: removal of map item means 
		 * that the file itself is kept. We start with a map of all files, then
		 * remove the map items that we want to keep in reality and finally 
		 * delete all remaining files.  
		 * 
		 * But there is a good reasons for it: its safe even when a new file is 
		 * added during this operation.   
		 */
		TreeMap<DateTime, File> filesToDelete = new TreeMap<>();
		// all files whose name can be parsed are potentially deleted
		parseFilenames(filesToDelete);
		// i.e. keep last n files
		removeLastItems(filesToDelete, ROTATION_COUNT);		
		// i.e. keep files of last n days		
		removeNewerThan(filesToDelete, ROTATION_DAYS);
		// i.e. keep first file of each month
		removeFirstOfEachMonth(filesToDelete);
		
		for (File file : filesToDelete.values()) {			
			logger.debug("backup rotation deletes file " + file);
			file.delete();
		}
		
		return filesToDelete.size();
	}


	/**
	 * Remove first map item (i.e. keep the file itself) of each month.
	 * 
	 * @param filesToDelete
	 */
	private void removeFirstOfEachMonth(TreeMap<DateTime, File> filesToDelete) {
		
		HashSet<DateTime> months = new HashSet<>();
		
		Iterator<DateTime> filesIter = filesToDelete.keySet().iterator();
		
		while (filesIter.hasNext()) {
			DateTime fileDate = filesIter.next();
			DateTime month = new DateTime(fileDate.getYear(), fileDate.getMonthOfYear(), 1, 0, 0);
			
			if (!months.contains(month)) {
				filesIter.remove();
				months.add(month);
			}
		}
	}
	
	/**
	 * Remove map items that are newer than the given date (i.e. keep the original files on disk).
	 * 
	 * @param filesToDelete
	 * @param count
	 */
	private void removeNewerThan(TreeMap<DateTime, File> filesToDelete, int days) {
		
		DateTime date = new DateTime();
		
		date = date.minusDays(days);
		
		Iterator<DateTime> filesIter = filesToDelete.keySet().iterator();
		
		while (filesIter.hasNext()) {
			DateTime fileDate = filesIter.next();
			
			if (fileDate.isAfter(date)) {
				filesIter.remove();
			}
		}
	}

	/**
	 * Remove last n files of the sorted map (i.e. keep the original files on disk).
	 * 
	 * @param filesToDelete
	 * @param count
	 */
	private void removeLastItems(TreeMap<DateTime, File> filesToDelete, int count) {
		
		for (int i = 0; i < count; i++) {
			if (filesToDelete.isEmpty()) {
				break;
			}
			filesToDelete.remove(filesToDelete.lastKey());
		}
	}

	/**
	 * Iterate over files in backupRoot and create a map of all files, whose filename starts with filenamePrefix,
	 * ends with filenamePostfix and has a legal DateTime string in between. The map is sorted in ascending order
	 * according the DateTime objects. 
	 * 
	 * @param filesToDelete
	 */
	private void parseFilenames(TreeMap<DateTime, File> filesToDelete) {
		
		for (File file : backupRoot.listFiles()) {
			String filename = file.getName();
			if (!filename.startsWith(filenamePrefix)) {
				continue;
			}
			
			if (!filename.endsWith(filenamePostfix)) {
				continue;
			}
			String dateString = filename.substring(filenamePrefix.length(), filename.length() - filenamePostfix.length());

			DateTime fileDate = null;
			try {
				fileDate = ISODateTimeFormat.dateTime().parseDateTime(dateString);			
			} catch (IllegalArgumentException e) {
				continue;
			}
			
			filesToDelete.put(fileDate, file);
		}
	}
	
	public File getBackupFile(DateTime date) {
		return new File(backupRoot, filenamePrefix + date + filenamePostfix);
	}
	
	public static void main(String[] args) throws IOException {
		
		File backupRoot = new File("test-backups");
		backupRoot.mkdir();
		
		BackupRotation rotation = new BackupRotation(backupRoot);
		
		rotation.createFakeBackupFiles();
						
		rotation.rotate();
		
		TreeMap<DateTime, File> keptFiles = new TreeMap<>();
		rotation.parseFilenames(keptFiles);
		
		for (File file : keptFiles.values()) {
			System.out.println(file);
		}
	}
	
	/**
	 * Create some empty files with backup-like filenames for testing.
	 * 
	 * @throws IOException
	 */
	private void createFakeBackupFiles() throws IOException {
		DateTime date = new DateTime(2010, 1, 1, 5, 00);
		
		while (date.isBeforeNow()) {
			
			File backupFile = getBackupFile(date);
			backupFile.createNewFile();
			
			date = date.plus(24*60*60*1000);						
		}
	}

}
