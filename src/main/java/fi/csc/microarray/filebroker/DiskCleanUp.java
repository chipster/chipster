package fi.csc.microarray.filebroker;

import java.io.File;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;

import fi.csc.microarray.util.Files;

public class DiskCleanUp {
	
	private static Logger logger = Logger.getLogger(DiskCleanUp.class);
	
	private File root;
	private int cleanUpTriggerLimitPercentage;
	private int cleanUpTargetPercentage;
	private int cleanUpMinimumFileAge;

	private final ExecutorService executor = Executors.newSingleThreadExecutor();
	private Future<?> lastCleanUp;
	private Object lastCleanUpLock = new Object(); // lock mustn't be null

	public DiskCleanUp(File root, int cleanUpTriggerLimitPercentage, int cleanUpTargetPercentage, int cleanUpMinimumFileAge) {
		this.root = root;
		this.cleanUpTriggerLimitPercentage = cleanUpTriggerLimitPercentage;
		this.cleanUpTargetPercentage = cleanUpTargetPercentage;
		this.cleanUpMinimumFileAge = cleanUpMinimumFileAge;
		
		logger.info("total space: " + FileUtils.byteCountToDisplaySize(root.getTotalSpace()));
		logger.info("usable space: " + FileUtils.byteCountToDisplaySize(root.getUsableSpace()));
		logger.info("cache clean up will start when usable space is less than: " + FileUtils.byteCountToDisplaySize((long) ((double)root.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100)) + " (" + (100-cleanUpTriggerLimitPercentage) + "%)");
		logger.info("cache clean target usable space is:  " + FileUtils.byteCountToDisplaySize((long) ((double)root.getTotalSpace()*(double)(100-cleanUpTargetPercentage)/100)) + " (" + (100-cleanUpTargetPercentage) + "%)");		
		logger.info("will not clean up files newer than: " + (cleanUpMinimumFileAge/3600) + "h");
	}
	
	public long getCleanUpTargetLimit() {
		return (long) ((double)root.getTotalSpace()*(double)(100-cleanUpTargetPercentage)/100);
	}

	public void scheduleCleanUp(final long requestedSize) {			
		try {
			runIfNotAlreadyRunning(new CleanUpRunnable(requestedSize), false);
		} catch (ExecutionException | InterruptedException e) {
			logger.warn("exception while cleaning cache", e);
		}		
	}
	
	public void cleanUpAndWait(long requestedSize) {
		try {
			runIfNotAlreadyRunning(new CleanUpRunnable(requestedSize), true);
		} catch (ExecutionException | InterruptedException e) {
			logger.warn("exception while cleaning cache", e);
		}
	}

	private void runIfNotAlreadyRunning(Runnable runnable, boolean wait) throws ExecutionException, InterruptedException {

		Future<?> cleanUp = null;
		
		synchronized (lastCleanUpLock) {			
			if (lastCleanUp == null || lastCleanUp.isDone()) {	
				lastCleanUp = executor.submit(runnable);
			} else {
				logger.info("cache cleanup already running, skipping this one");
			}
			cleanUp = lastCleanUp;
		}
		
		if (wait) {
			cleanUp.get();
		}
	}

	public class CleanUpRunnable implements Runnable {
		private long requestedSize;

		public CleanUpRunnable(long requestedSize) {
			this.requestedSize = requestedSize;
		}
		public void run () {

			long cleanUpBeginTime = System.currentTimeMillis();
			long cleanUpTargetLimit = getCleanUpTargetLimit();
			logger.info("cache cleanup, target usable space: " + FileUtils.byteCountToDisplaySize(requestedSize + cleanUpTargetLimit) + 
					" (" + FileUtils.byteCountToDisplaySize(requestedSize) + " + " + FileUtils.byteCountToDisplaySize(cleanUpTargetLimit) + 
					" (" + (100-cleanUpTargetPercentage) + "%)");
			Files.makeSpaceInDirectory(root, requestedSize + cleanUpTargetLimit, cleanUpMinimumFileAge, TimeUnit.SECONDS);
			logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms, usable space now " + FileUtils.byteCountToDisplaySize(root.getUsableSpace())); 
		}
	}

	/**
	 * @return requestedSize is available
	 */
	public boolean scheduleCleanUpIfNecessary(long requestedSize) {
		long usableSpaceSoftLimit =  (long) ((double)root.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100);		
		logger.debug("usable space soft limit is: " + usableSpaceSoftLimit);
			
		boolean spaceAvailable;
		
		if (root.getUsableSpace() - requestedSize >= usableSpaceSoftLimit) {
			// space available, clean up limit will not be reached
			logger.debug("enough space available, no need to do anything");
			spaceAvailable = true;
			
		} else {
			// space available, clean up soft limit will be reached
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(requestedSize) + " usable: " + FileUtils.byteCountToDisplaySize(root.getUsableSpace()) + 
					", usable space soft limit: " + FileUtils.byteCountToDisplaySize(usableSpaceSoftLimit) + " (" + (100-cleanUpTriggerLimitPercentage) + 
					"%) will be reached --> scheduling clean up");
			spaceAvailable = true;
			
			// schedule clean up
			scheduleCleanUp(requestedSize);
		}
		return spaceAvailable;
	}
}
