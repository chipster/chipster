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
	private long minimumSpaceForAcceptUpload;

	private final ExecutorService executor = Executors.newSingleThreadExecutor();
	private Future<?> lastCleanUp;
	private Object lastCleanUpLock = new Object(); // lock mustn't be null


	public DiskCleanUp(File root, int cleanUpTriggerLimitPercentage, int cleanUpTargetPercentage, int cleanUpMinimumFileAge, long minimumSpaceForAcceptUpload) {
		this.root = root;
		this.cleanUpTriggerLimitPercentage = cleanUpTriggerLimitPercentage;
		this.cleanUpTargetPercentage = cleanUpTargetPercentage;
		this.cleanUpMinimumFileAge = cleanUpMinimumFileAge;
		this.minimumSpaceForAcceptUpload = minimumSpaceForAcceptUpload;
		
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
			
			//FIXME
			try {
				Thread.sleep(5*60_000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			//FIXME

			long cleanUpBeginTime = System.currentTimeMillis();
			long cleanUpTargetLimit = getCleanUpTargetLimit();
			logger.info("cache cleanup, target usable space: " + FileUtils.byteCountToDisplaySize(requestedSize + cleanUpTargetLimit) + 
					" (" + FileUtils.byteCountToDisplaySize(requestedSize) + " + " + FileUtils.byteCountToDisplaySize(cleanUpTargetLimit) + 
					", " + (100-cleanUpTargetPercentage) + "%)");
			Files.makeSpaceInDirectory(root, requestedSize + cleanUpTargetLimit, cleanUpMinimumFileAge, TimeUnit.SECONDS);
			logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms, usable space now " + FileUtils.byteCountToDisplaySize(root.getUsableSpace())); 
		}
	}

	/**
	 * Handle space request.
	 * 
	 * This will start a clean up when necessary. On a soft limit, the current
	 * request can be satisfied immediately and the clean up is only scheduled.
	 * If a hard limit is reached and allowWait is true, this will wait until
	 * the clean up is done.
	 * 
	 * @param requestedSize
	 * @param allowWait
	 *            When the space request can't be satisfied immediately, should
	 *            we wait for clean up or just return false.
	 * @return true if requested bytes are available
	 */
	public boolean spaceRequest(long requestedSize, boolean allowWait) {
		long usableSpaceSoftLimit =  (long) ((double)root.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100);		
		long usableSpaceHardLimit = minimumSpaceForAcceptUpload;	
		
		// deal with the weird config case of soft limit being smaller than hard limit
		if (usableSpaceSoftLimit < usableSpaceHardLimit) {
			usableSpaceSoftLimit = usableSpaceHardLimit;
		}
		
		logger.debug("usable space soft limit is: " + usableSpaceSoftLimit);
		logger.debug("usable space hard limit is: " + usableSpaceHardLimit);
		
		boolean spaceAvailable;
			
		long availableAfterCleanUp = getCleanUpTargetLimit() - minimumSpaceForAcceptUpload;
		
		if (getBytesAvailable() > usableSpaceSoftLimit) {
			
			// space available, clean up limit will not be reached
			logger.debug("enough space available, no need to do anything");
			spaceAvailable = true;
			
		} else if (getBytesAvailable() > requestedSize) {

			// space available, but clean up soft limit will be reached
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(requestedSize) + " usable: " + FileUtils.byteCountToDisplaySize(root.getUsableSpace()) + 
					", usable space soft limit: " + FileUtils.byteCountToDisplaySize(usableSpaceSoftLimit) + " (" + (100-cleanUpTriggerLimitPercentage) + 
					"%) will be reached --> scheduling clean up");

			// schedule clean up
			scheduleCleanUp(requestedSize);
			spaceAvailable = true;

		} else if (availableAfterCleanUp > requestedSize){

			// will run out of usable space, wait for clean up		
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(requestedSize) + " usable: " + FileUtils.byteCountToDisplaySize(root.getUsableSpace()) + 
					", not enough space --> wait for clean up");

			if (allowWait) {
				cleanUpAndWait(requestedSize);
			} else {
				scheduleCleanUp(requestedSize);
				logger.info("space request can't be satisfied, but waiting isn't allowed");
			}

			logger.info("not accepting upload if less than " + FileUtils.byteCountToDisplaySize(minimumSpaceForAcceptUpload) + " usable space after upload");

			// check if cleaned up enough 
			if (getBytesAvailable() > requestedSize) {
				logger.info("enough space after cleaning");
				spaceAvailable = true;
			} else {
				logger.info("not enough space after cleaning");
				spaceAvailable = false;
			}			

		} else {
			// request more than total, no can do
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(requestedSize) + ", usable: " + FileUtils.byteCountToDisplaySize(root.getUsableSpace()) + 
					", maximum space: " + FileUtils.byteCountToDisplaySize(root.getTotalSpace()) + 
					", minimum usable: " + FileUtils.byteCountToDisplaySize(minimumSpaceForAcceptUpload) + 
					" --> not possible to make enough space");

			spaceAvailable = false;
		}
		
		return spaceAvailable;
	}
	
	private long getBytesAvailable() {
		return root.getUsableSpace() - minimumSpaceForAcceptUpload;
	}
}
