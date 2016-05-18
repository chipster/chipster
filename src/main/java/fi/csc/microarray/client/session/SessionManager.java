package fi.csc.microarray.client.session;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.UUID;

import javax.swing.Timer;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.filebroker.DerbyMetadataServer;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.filebroker.QuotaExceededException;
import fi.csc.microarray.messaging.AuthCancelledException;
import fi.csc.microarray.messaging.admin.StorageAdminAPI.StorageEntryMessageListener;
import fi.csc.microarray.security.CryptoKey;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.Files;

public class SessionManager {

	private Logger logger = Logger.getLogger(SessionManager.class);

	/**
	 * Reports session validation related problems.
	 */
	@SuppressWarnings("serial")
	public static class ValidationException extends Exception {
		public ValidationException(String validationDetails) {
			super(validationDetails);
		}
	}

	public interface SessionManagerCallback {
		public void showDialog(String title, String message, String details,
				Severity severity, boolean modal,
				DetailsVisibility detailsVisibility);

		public void reportException(Exception e);

		public void sessionChanged(SessionChangedEvent e);

		public List<OperationRecord> getUnfinishedJobs();

		public void continueJobs(List<OperationRecord> unifinishedJobs);
	}

	protected static final String ALIVE_SIGNAL_FILENAME = "i_am_alive";
	protected static final int SESSION_BACKUP_INTERVAL = 5 * 1000;

	private String sessionNotes;
	private String currentSessionName;
	private String currentRemoteSession;

	protected boolean unsavedChanges = false;
	protected boolean unbackuppedChanges = false;

	protected File aliveSignalFile;
	private LinkedList<File> deadDirectories = new LinkedList<File>();

	private DataManager dataManager;
	private FileBrokerClient fileBrokerClient;
	private SessionManagerCallback callback;
	private TaskExecutor taskExecutor;

	/**
	 * @param dataManager
	 * @param taskExecutor 
	 * @param fileBrokerClient
	 * @param callback
	 *            if null, SessionChangedEvents are ignored and error messages
	 *            are thrown as RuntimeExceptions
	 * @throws IOException
	 */
	public SessionManager(final DataManager dataManager,
			TaskExecutor taskExecutor, FileBrokerClient fileBrokerClient, SessionManagerCallback callback)
			throws IOException {
		this.dataManager = dataManager;
		this.taskExecutor = taskExecutor;
		this.fileBrokerClient = fileBrokerClient;
		if (callback == null) {
			this.callback = new BasicSessionManagerCallback();
		} else {
			this.callback = callback;
		}

		// Remember changes to confirm close only when necessary and to backup
		// when necessary
		dataManager.addDataChangeListener(new DataChangeListener() {
			public void dataChanged(DataChangeEvent event) {
				unsavedChanges = true;
				unbackuppedChanges = true;
			}
		});

		// Start checking if background backup is needed
		aliveSignalFile = new File(dataManager.getRepository(), "i_am_alive");
		aliveSignalFile.createNewFile();
		aliveSignalFile.deleteOnExit();

		Timer timer = new Timer(SESSION_BACKUP_INTERVAL, new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				aliveSignalFile.setLastModified(System.currentTimeMillis()); // touch
																				// the
																				// file
				if (unbackuppedChanges) {

					File sessionFile = UserSession.findBackupFile(
							dataManager.getRepository(), true);
					sessionFile.deleteOnExit();

					try {
						saveLightweightSession(sessionFile);

					} catch (Exception e1) {
						logger.warn(e1); // do not care that much about failing
											// session backups
					}
				}
				unbackuppedChanges = false;
			}
		});

		timer.setCoalesce(true);
		timer.setRepeats(true);
		timer.setInitialDelay(SESSION_BACKUP_INTERVAL);
		timer.start();
	}

	public class BasicSessionManagerCallback implements SessionManagerCallback {

		@Override
		public void showDialog(String title, String message, String details,
				Severity severity, boolean modal,
				DetailsVisibility detailsVisibility) {
			throw new RuntimeException(title + ": " + message + " (" + details
					+ ")");
		}

		@Override
		public void reportException(Exception e) {
			throw new RuntimeException(e);
		}

		@Override
		public void sessionChanged(SessionChangedEvent e) {
			// ignore
		}

		@Override
		public List<OperationRecord> getUnfinishedJobs() {
			return new ArrayList<>();
		}

		@Override
		public void continueJobs(List<OperationRecord> unfinishedJobs) {
			logger.info("Continuation of jobs isn't implemented in "
					+ this.getClass().getSimpleName());
		}
	}

	public void loadLocalSession(File sessionFile, boolean isDataless)
			throws Exception {
		loadLocalSession(sessionFile, isDataless, null);
	}

	/**
	 * Load session from a file.
	 * 
	 * @param xOffset
	 * 
	 * @see #saveSession(File, ClientApplication)
	 */
	public void loadLocalSession(File sessionFile, boolean isDataless,
			Integer xOffset) throws Exception {
		SessionLoader sessionLoader = new SessionLoader(sessionFile,
				isDataless, dataManager);
		loadSession(sessionLoader, xOffset);
	}

	public void loadStorageSession(String sessionId) throws Exception {
		loadStorageSession(sessionId, null);
	}

	/**
	 * Load remote session from an URL.
	 * 
	 * @param xOffset
	 * 
	 * @see #saveStorageSession(String)
	 */
	public void loadStorageSession(String sessionId, Integer xOffset)
			throws Exception {
		SessionLoader sessionLoader = new SessionLoader(sessionId, dataManager);
		loadSession(sessionLoader, xOffset);
	}

	private void loadSession(SessionLoader sessionLoader, Integer xOffset)
			throws Exception {
		sessionLoader.setXOffset(xOffset);
		sessionLoader.loadSession();
		// set session notes directly to field because this doesn't trigger
		// unsaved changes
		sessionNotes = sessionLoader.getSessionNotes();
		callback.continueJobs(sessionLoader.getUnifinishedJobs());
	}

	/**
	 * Saves session (all data: beans, folder structure, operation metadata,
	 * links etc.) to a file. File is a zip file with all the data files and one
	 * metadata file.
	 * 
	 * @param unfinishedJobs
	 * 
	 * @return true if the session was saved perfectly
	 * @throws Exception
	 */
	public void saveSession(File sessionFile) throws Exception {

		// save session file
		boolean metadataValid = false;
		SessionSaver sessionSaver = new SessionSaver(sessionFile, dataManager);
		sessionSaver.setSessionNotes(sessionNotes);
		sessionSaver.setUnfinishedJobs(callback.getUnfinishedJobs());
		metadataValid = sessionSaver.saveSession();

		// check validation
		if (!metadataValid) {
			// save was successful but metadata validation failed, file might be
			// usable
			String validationDetails = sessionSaver.getValidationErrors();
			throw new ValidationException(validationDetails);
		}
	}

	/**
	 * Saves lightweight session (folder structure, operation metadata, links
	 * etc.) to a file. Does not save actual data inside databeans.
	 * 
	 * @param unfinishedJobs
	 * 
	 * @return true if the session was saved perfectly
	 * @throws Exception
	 */
	public void saveLightweightSession(File sessionFile) throws Exception {

		SessionSaver sessionSaver = new SessionSaver(sessionFile, dataManager);
		sessionSaver.setSessionNotes(sessionNotes);
		sessionSaver.setUnfinishedJobs(callback.getUnfinishedJobs());
		sessionSaver.saveLightweightSession();
	}

	/**
	 * Returns debug print out of current session state.
	 * 
	 * @return print out of session state
	 */
	public String printSession() {
		StringBuffer buffer = new StringBuffer();
		SessionSaver.dumpSession(dataManager.getRootFolder(), buffer);
		return buffer.toString();
	}

	public String saveStorageSession(String name) throws Exception {

		String sessionId = CryptoKey.generateRandom();
		SessionSaver sessionSaver = new SessionSaver(sessionId, dataManager);
		sessionSaver.setSessionNotes(sessionNotes);
		sessionSaver.setUnfinishedJobs(callback.getUnfinishedJobs());
		// upload/move data files and upload metadata files, if needed
		LinkedList<String> dataIds = sessionSaver.saveStorageSession();

		// add metadata to file broker database (make session visible)
		fileBrokerClient.saveRemoteSession(name, sessionId, dataIds);

		return sessionId;
	}

	/**
	 * Like saveStorageSession, but upload necessary files to cache instead of
	 * storage and returns sessionId instead of writing it to metadata database.
	 * 
	 * @param unfinishedJobs
	 * 
	 * @return
	 * 
	 * @return sessionId
	 * @throws Exception
	 */
	public String saveFeedbackSession() throws Exception {

		String sessionId = CryptoKey.generateRandom();
		SessionSaver sessionSaver = new SessionSaver(sessionId, dataManager);
		sessionSaver.setSessionNotes(sessionNotes);
		sessionSaver.setUnfinishedJobs(callback.getUnfinishedJobs());
		// upload/move data files and upload metadata files, if needed
		sessionSaver.saveFeedbackSession();

		return sessionId;
	}

	public List<DbSession> listRemoteSessions() throws FileBrokerException, AuthCancelledException {
		return fileBrokerClient.listRemoteSessions();
	}
	
	public StorageEntryMessageListener getStorageUsage() throws FileBrokerException, InterruptedException {
		return fileBrokerClient.getStorageUsage();
	}

	public void setSession(File sessionFile, String sessionId)
			throws MalformedURLException, FileBrokerException, AuthCancelledException {
		if (sessionFile != null) {
			currentRemoteSession = null;
			String oldValue = currentSessionName;
			currentSessionName = sessionFile.getName().replace(".zip", "");
			callback.sessionChanged(new SessionChangedEvent(this, "session",
					oldValue, currentSessionName));
		} else if (sessionId != null) {
			String oldValue = currentRemoteSession;
			currentRemoteSession = sessionId;
			DbSession session = findSessionWithUuid(listRemoteSessions(), sessionId);
			if (session != null) {
				currentSessionName = session.getName();
			} else {
				// help request sessions don't have a name
				currentSessionName = null;
			}
			callback.sessionChanged(new SessionChangedEvent(this, "session",
					oldValue, currentRemoteSession));
		} else {
			String oldValue = currentRemoteSession != null ? currentRemoteSession
					: currentSessionName;
			currentRemoteSession = null;
			currentSessionName = null;
			callback.sessionChanged(new SessionChangedEvent(this, "session",
					oldValue, null));
		}
	}

	public String getSessionName() {
		return currentSessionName;
	}

	public DbSession findSessionWithName(List<DbSession> sessions, String name) throws MalformedURLException {
		for (DbSession session : sessions) {
			if (session.getName() != null && session.getName().equals(name)) {
				return session;
			}
		}
		return null;
	}
	
	public DbSession findSessionWithUuid(List<DbSession> sessions, String uuid) throws MalformedURLException {
		for (DbSession session : sessions) {
			if (session.getDataId() != null && session.getDataId().equals(uuid)) {
				return session;
			}
		}
		return null;
	}

	public String getSessionNotes() {
		return this.sessionNotes;
	}

	public void setSessionNotes(String notes) {
		if (this.sessionNotes != null && !this.sessionNotes.equals(notes)) {
			unsavedChanges = true;
		}
		this.sessionNotes = notes;
	}

	public boolean isCurrentRemoteSession(String sessionUuid) {
		return currentRemoteSession != null
				&& currentRemoteSession.equals(sessionUuid);
	}

	public static class SessionChangedEvent extends PropertyChangeEvent {

		public SessionChangedEvent(Object source, String propertyName,
				Object oldValue, Object newValue) {
			super(source, propertyName, oldValue, newValue);
		}
	}

	public boolean areCloudSessionsEnabled() {
		boolean conf = DirectoryLayout.getInstance().getConfiguration()
				.getBoolean("client", "enable-cloud-sessions");
		boolean specialUser = DerbyMetadataServer.DEFAULT_EXAMPLE_SESSION_OWNER
				.equals(Session.getSession().getUsername());

		return conf || specialUser;
	}

	public void restoreSessionAndWait(File file) {
		loadSessionAndWait(file, null, true, true, false, 0);
	}

	public void loadSessionAndWait(final File sessionFile,
			final String sessionId, final boolean isDataless,
			final boolean clearDeadTempDirs, final boolean isExampleSession) {
		loadSessionAndWait(sessionFile, sessionId, isDataless,
				clearDeadTempDirs, isExampleSession, null);
	}

	public void loadSessionAndWait(final File sessionFile,
			final String sessionId, final boolean isDataless,
			final boolean clearDeadTempDirs, final boolean isExampleSession,
			Integer xOffset) {

		// check that it's a valid session file
		if (!isDataless) {
			if (!UserSession.isValidSessionFile(sessionFile)) {
				Session.getSession()
						.getApplication()
						.showDialog("Could not open session file.",
								"The given file is not a valid session file.",
								"", Severity.INFO, true);
				return;
			}
		}

		/*
		 * If there wasn't data or it was just cleared, there is no need to warn
		 * about saving after opening session. However, if there was datasets
		 * already, combination of them and new session can be necessary to
		 * save. This has to set after the import.
		 */
		boolean somethingToSave = dataManager.databeans().size() != 0;

		try {
			if (sessionFile != null) {
				loadLocalSession(sessionFile, isDataless);
			} else {
				loadStorageSession(sessionId);
			}

			setSession(sessionFile, sessionId);

		} catch (Exception e) {
			if (isExampleSession) {
				Session.getSession()
						.getApplication()
						.showDialog(
								"Opening example session failed.",
								"Please restart "
										+ Session.getSession()
												.getPrimaryModule()
												.getDisplayName()
										+ " to update example session links or see the details for more information.",
								Exceptions.getStackTrace(e), Severity.INFO,
								true, DetailsVisibility.DETAILS_HIDDEN, null);
			} else {
				Session.getSession()
						.getApplication()
						.showDialog(
								"Opening session failed.",
								"Unfortunately the session could not be opened properly. Please see the details for more information.",
								Exceptions.getStackTrace(e), Severity.WARNING,
								true, DetailsVisibility.DETAILS_HIDDEN, null);
			}
			logger.error("loading session failed", e);
		}

		unsavedChanges = somethingToSave;

		// If this was restored session, clear dead temp directories in the end.
		// It is done inside this method to avoid building synchronization
		// between
		// session loading and temp directory cleaning during restore.
		if (clearDeadTempDirs) {
			clearDeadTempDirectories();
		}
	}

	public boolean saveSessionAndWait(boolean isRemote, File localFile,
			String remoteSessionName) {

		try {
			String sessionId = null;

			if (isRemote) {
				sessionId = saveStorageSession(remoteSessionName);
			} else {
				saveSession(localFile);
			}

			setSession(localFile, sessionId);

			unsavedChanges = false;
			return true;

		} catch (ValidationException e) {
			callback.showDialog(
					"Problem with saving the session",
					"All the datasets were saved successfully, but there were troubles with saving "
							+ "the session information about them. This means that there may be problems when "
							+ "trying to open the saved session file later on.\n"
							+ "\n"
							+ "If you have important unsaved "
							+ "datasets in this session, it might be a good idea to export such datasets using the "
							+ "File -> Export functionality.", e.getMessage(),
					Severity.WARNING, true, DetailsVisibility.DETAILS_HIDDEN);

			return false;

		} catch (QuotaExceededException e) {
			callback.showDialog(
					"Quota exceeded",
					"Saving session failed, because your disk space quota was exceeded.\n"
							+ "\n"
							+ "Please contact server maintainers to apply for more quota, remove some old sessions "
							+ "to free more disk space or save the session on your computer using the "
							+ "File -> Save local session functionality. ",
					e.getMessage(), Severity.WARNING, true,
					DetailsVisibility.DETAILS_ALWAYS_HIDDEN);
			return false;

		} catch (Exception e) {
			callback.showDialog(
					"Saving session failed",
					"Unfortunately your session could not be saved. Please see the details for more "
							+ "information.\n"
							+ "\n"
							+ "If you have important unsaved datasets in this session, it might be "
							+ "a good idea to export such datasets using the File -> Export functionality.",
					Exceptions.getStackTrace(e), Severity.WARNING, true,
					DetailsVisibility.DETAILS_HIDDEN);
			return false;
		}
	}

	public void clearSessionWithoutConfirming() throws MalformedURLException,
	FileBrokerException, AuthCancelledException {
		dataManager.deleteAllDataItems();
		taskExecutor.clear();
		setSessionNotes(null);
		setSession(null, null);
		unsavedChanges = false;
	}

	public boolean removeRemoteSession(String sessionUuid) throws FileBrokerException, AuthCancelledException {

		if (currentRemoteSession != null
				&& currentRemoteSession.equals(sessionUuid)
				&& !dataManager.databeans().isEmpty()) {
			callback.showDialog(
					"Remove prevented",
					"You were trying to remove a cloud session that is your last saved session. "
							+ "Removal of this session is prevented, because it may be the only copy of your current "
							+ "datasets. If you want to keep the datasets, please save them as a sessions first. If you want to remove "
							+ "the datasets, please delete them before removing the cloud session.",
					null, Severity.INFO, true,
					DetailsVisibility.DETAILS_ALWAYS_HIDDEN);
			return false;
		}

		fileBrokerClient.removeRemoteSession(sessionUuid);
		return true;
	}

	public boolean hasUnsavedChanges() {
		return unsavedChanges;
	}

	public void clearDeadTempDirectories() {

		// Try to clear dead temp directories
		try {
			for (File dir : deadDirectories) {
				Files.delTree(dir);
			}
		} catch (Exception e) {
			callback.reportException(e);
		}

		// Remove them from bookkeeping in any case
		deadDirectories.clear();
	}

	/**
	 * Collects all dead temp directories and returns the most recent that has a
	 * restorable session .
	 */
	public File checkTempDirectories() throws IOException {

		Iterable<File> tmpDirectories = dataManager.listAllRepositories();
		File mostRecentDeadSignalFile = null;

		for (File directory : tmpDirectories) {

			// Skip current temp directory
			if (directory.equals(dataManager.getRepository())) {
				continue;
			}

			// Check is it alive, wait until alive file should have been updated
			File aliveSignalFile = new File(directory, ALIVE_SIGNAL_FILENAME);
			long originalLastModified = aliveSignalFile.lastModified();
			boolean unsuitable = false;
			while ((System.currentTimeMillis() - aliveSignalFile.lastModified()) < 2 * SESSION_BACKUP_INTERVAL) {

				// Updated less than twice the interval time ago
				// ("not too long ago"), so keep on checking
				// until we see new update that confirms it is alive, or have
				// waited long
				// enough that the time since last update grows larger than
				// twice the interval.

				// Check if restorable
				if (UserSession.findBackupFile(directory, false) == null) {
					// Does not have backup file, so not interesting for backup.
					// Should be removed anyway, but removing empty directories
					// is not
					// important enough to warrant the extra waiting that
					// follows next.
					// So we will skip this and if it was dead, it will be
					// anyway
					// cleaned away in the next client startup.

					unsuitable = true;
					break;
				}

				// Check if updated
				if (aliveSignalFile.lastModified() != originalLastModified) {
					unsuitable = true;
					break; // we saw an update, it is alive
				}

				// Wait for it to update
				try {
					Thread.sleep(1000); // 1 second
				} catch (InterruptedException e) {
					// ignore
				}
			}

			if (!unsuitable) {
				// It is dead, might be the one that should be recovered, check
				// that
				deadDirectories.add(directory);
				File deadSignalFile = new File(directory, ALIVE_SIGNAL_FILENAME);
				if (UserSession.findBackupFile(directory, false) != null
						&& (mostRecentDeadSignalFile == null || mostRecentDeadSignalFile
								.lastModified() < deadSignalFile.lastModified())) {

					mostRecentDeadSignalFile = deadSignalFile;

				}
			}
		}

		return mostRecentDeadSignalFile != null ? mostRecentDeadSignalFile
				.getParentFile() : null;
	}

	public void setUnsavedChanges() {
		this.unsavedChanges = true;
	}
	
	public UUID getSessionId() {
		if (currentRemoteSession != null) {
			return UUID.fromString(currentRemoteSession);
		} else {
			return null;
		}
	}
}
