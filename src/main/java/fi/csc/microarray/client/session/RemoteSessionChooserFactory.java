package fi.csc.microarray.client.session;

import java.awt.Dimension;
import java.io.File;
import java.net.MalformedURLException;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;
import javax.swing.JFileChooser;

import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.serverfiles.ServerFile;
import fi.csc.microarray.client.serverfiles.ServerFileSystemView;
import fi.csc.microarray.client.serverfiles.ServerFileUtils;
import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.filebroker.DerbyMetadataServer;
import fi.csc.microarray.filebroker.FileBrokerException;

public class RemoteSessionChooserFactory {
	
	private SwingClientApplication app;

	public RemoteSessionChooserFactory(SwingClientApplication app) {
		this.app = app;
	}

	/**
	 * Excludes example sessions.
	 * 
	 * @return
	 * @throws JMSException
	 * @throws Exception
	 * @throws MalformedURLException
	 */
	private JFileChooser populateFileChooserFromServer() throws JMSException, Exception, MalformedURLException {
		return populateFileChooserFromServer(true);
	}
	
	
	private JFileChooser populateFileChooserFromServer(boolean excludeExampleSessions) throws JMSException, Exception, MalformedURLException {
		JFileChooser sessionFileChooser = new JFileChooser();
		sessionFileChooser.setMultiSelectionEnabled(false);
		updateRemoteSessions(app.getSessionManager(), sessionFileChooser, excludeExampleSessions);
		app.fixFileChooserFontSize(sessionFileChooser);
		return sessionFileChooser;
	}

	

	public static ServerFileSystemView updateRemoteSessions(SessionManager sessionManager,
			JFileChooser sessionFileChooser, boolean excludeExampleSessions) throws FileBrokerException, MalformedURLException {
		List<DbSession> sessions = sessionManager.listRemoteSessions();
		
		// hide example sessions
		if (excludeExampleSessions) {
			List<DbSession> sessionsToHide = new LinkedList<DbSession>();
			for (DbSession session: sessions) {
				if (session.getName().startsWith(DerbyMetadataServer.DEFAULT_EXAMPLE_SESSION_FOLDER)) {
					sessionsToHide.add(session);
				}
			}
			sessions.removeAll(sessionsToHide);
		}
			
		ServerFileSystemView view = ServerFileSystemView.parseFromPaths(ServerFile.SERVER_SESSION_ROOT_FOLDER, sessions);		
		sessionFileChooser.setFileSystemView(view);
		// without this the GUI doesn't update when a session is removed
		sessionFileChooser.setCurrentDirectory(null);
		sessionFileChooser.setCurrentDirectory(view.getRoot());		
		sessionFileChooser.putClientProperty("sessions", sessions);
		// without this removed sessions are shown after folder change
		sessionFileChooser.updateUI();
		return view;
	}

	public JFileChooser getExampleSessionChooser() throws MalformedURLException, JMSException, Exception {
		JFileChooser exampleSessionFileChooser = populateFileChooserFromServer(false);
		exampleSessionFileChooser.setSelectedFile(new File("session"));
		ServerFileUtils.hideJFileChooserButtons(exampleSessionFileChooser);
		exampleSessionFileChooser.setPreferredSize(new Dimension(800, 600));

		ServerFileSystemView view = (ServerFileSystemView) exampleSessionFileChooser.getFileSystemView();
		exampleSessionFileChooser.setCurrentDirectory(view.getExampleSessionDir());				
		
		return exampleSessionFileChooser;
	}

	public JFileChooser getRemoteSessionChooser() throws MalformedURLException, JMSException, Exception {
		JFileChooser remoteSessionFileChooser = populateFileChooserFromServer();
		remoteSessionFileChooser.setSelectedFile(new File("session"));
		remoteSessionFileChooser.setPreferredSize(new Dimension(800, 600));
		remoteSessionFileChooser.setAccessory(new RemoteSessionAccessory(remoteSessionFileChooser, app.getSessionManager(), app));
		ServerFileUtils.hideJFileChooserButtons(remoteSessionFileChooser);
			
		return remoteSessionFileChooser;
	}

	public JFileChooser getManagementChooser() {

		JFileChooser sessionFileChooser = null;

		try {
			// fetch current sessions to show in the dialog and create it
			sessionFileChooser = populateFileChooserFromServer();

		} catch (Exception e) {
			throw new RuntimeException(e);
		}


		// tune GUI
		sessionFileChooser.setDialogTitle("Manage");

		sessionFileChooser.setPreferredSize(new Dimension(800, 600));
		sessionFileChooser.setAccessory(new RemoteSessionAccessory(sessionFileChooser, app.getSessionManager(), app));		

		// hide buttons that we don't need
		ServerFileUtils.hideJFileChooserButtons(sessionFileChooser);
		ServerFileUtils.hideApproveButton(sessionFileChooser);
		ServerFileUtils.setCancelButtonText(sessionFileChooser, "Close");

		return sessionFileChooser;
	}
}
