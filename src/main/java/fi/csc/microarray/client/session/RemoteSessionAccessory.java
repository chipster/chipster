package fi.csc.microarray.client.session;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.net.MalformedURLException;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.jms.JMSException;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JProgressBar;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.border.LineBorder;

import net.miginfocom.swing.MigLayout;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.operation.ColoredCircleIcon;
import fi.csc.microarray.client.serverfiles.ServerFile;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.filebroker.DerbyMetadataServer;
import fi.csc.microarray.messaging.admin.StorageAdminAPI.StorageEntryMessageListener;
import fi.csc.microarray.messaging.admin.StorageEntry;
import fi.csc.microarray.util.Strings;

public class RemoteSessionAccessory extends JPanel implements ActionListener, PropertyChangeListener {

	private static final Color LOW_DISK_USAGE_COLOR = VisualConstants.COLOR_BLUE_GREEN;
	private static final Color HIGH_DISK_USAGE_COLOR = VisualConstants.COLOR_ORANGE;
	
	private String disclaimer = ""
			+ "Storage here is for working copies and may not be backed up. Store "
			+ "another copy of all your valuable data elsewhere.";

	private JPanel panel = new JPanel();
	private JLabel manageTitle = new JLabel("Selected session");
	private JLabel diskUsageTitle = new JLabel("Disk usage");
	private JLabel disclaimerTitle = new JLabel("No backups");
	private JLabel previewLabel = new JLabel(" ");
	private JButton removeButton = new JButton("Remove");
	private JProgressBar quotaBar = new JProgressBar();
	private JTextArea disclaimerText = new JTextArea(disclaimer);
	private JLabel lowDiskUsageIcon = new JLabel(new ColoredCircleIcon(LOW_DISK_USAGE_COLOR));
	private JLabel highDiskUsageIcon = new JLabel(new ColoredCircleIcon(HIGH_DISK_USAGE_COLOR));	
	private JTextArea lowDiskUsageText = new JTextArea();
	private JTextArea highDiskUsageText = new JTextArea();
	private JFileChooser fileChooser;
	private SessionManager sessionManager;
	private SwingClientApplication app;
	
	private ExecutorService previewExecutor = Executors.newFixedThreadPool(1);

	public RemoteSessionAccessory(JFileChooser fileChooser, SessionManager sessionManager, SwingClientApplication app) {
		
		this.fileChooser = fileChooser;
		this.sessionManager = sessionManager;
		this.app = app;
		
		fileChooser.addPropertyChangeListener(this);		
		
		panel.setLayout(new MigLayout("", "[fill]", ""));
		panel.setBackground(Color.white);
		panel.setBorder(new LineBorder(VisualConstants.TOOL_LIST_BORDER_COLOR));		

		removeButton.addActionListener(this);
		
		disclaimerText.setLineWrap(true);
		lowDiskUsageText.setLineWrap(true);
		highDiskUsageText.setLineWrap(true);
		
		disclaimerText.setWrapStyleWord(true);
		lowDiskUsageText.setWrapStyleWord(true);
		highDiskUsageText.setWrapStyleWord(true);
		
		disclaimerText.setEditable(false);
		lowDiskUsageText.setEditable(false);
		highDiskUsageText.setEditable(false);
		
		disclaimerText.setOpaque(false);
		lowDiskUsageText.setOpaque(false);
		highDiskUsageText.setOpaque(false);
		
		quotaBar.setStringPainted(true);
		quotaBar.setBackground(Color.white);
		// get rid of a blue shadow by removing the border of the JProggresBar
		// and using JPanel to show the border
		quotaBar.setBorderPainted(false);
		quotaBar.setBorder(null);
		JPanel quotaPanel = new JPanel(new MigLayout("fill, gap 0!, insets 0"));
		quotaPanel.setBorder(new LineBorder(Color.GRAY));
		quotaPanel.add(quotaBar, "growx, width 300px");
		
		manageTitle.setFont(UIManager.getFont("TitledBorder.font"));
		diskUsageTitle.setFont(UIManager.getFont("TitledBorder.font"));
		disclaimerTitle.setFont(UIManager.getFont("TitledBorder.font"));
		
		manageTitle.setForeground(UIManager.getColor("TitledBorder.titleColor"));
		diskUsageTitle.setForeground(UIManager.getColor("TitledBorder.titleColor"));		
		disclaimerTitle.setForeground(UIManager.getColor("TitledBorder.titleColor"));

		panel.add(manageTitle, "span, wrap");
		panel.add(previewLabel, "skip, wrap");
		panel.add(removeButton, "sizegroup actions, skip, growx 0, wrap");
		
		panel.add(diskUsageTitle, "span, wrap");
		panel.add(quotaPanel, "skip, wrap");
		panel.add(lowDiskUsageIcon, "aligny top, gap top 5");
		panel.add(lowDiskUsageText, "wrap");
		panel.add(highDiskUsageIcon, "aligny top, gap top 5");
		panel.add(highDiskUsageText, "wrap");
				
		panel.add(disclaimerTitle, "span, wrap");
		panel.add(disclaimerText, "skip, wrap");
		
		this.setLayout(new MigLayout("fill, insets 0 0 1 1"));
		this.add(panel, "grow");

		update();
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		if (e.getSource() == removeButton) {
			removeSelectedSession();
		} 
	}

	private void removeSelectedSession() {
		
		String sessionUuid = null;
		
		try {
			sessionUuid = getSelectedSession().getDataId();
			
			if (sessionUuid == null) {
				throw new RuntimeException("session not found");
			}
		} catch (Exception e) {
			throw new RuntimeException("internal error: URL or name from save dialog was invalid"); // should never happen
		}
		
		try {
			// remove selected session
			if (sessionManager.removeRemoteSession(sessionUuid)) {			
				
				update();
			}

		} catch (JMSException e) {
			app.reportException(e);
		}
	}

	/**
	 * @return UUID of the selected session or null if no session is selected
	 * @throws MalformedURLException
	 */
	private DbSession getSelectedSession() throws MalformedURLException {
		File selectedFile = fileChooser.getSelectedFile();
		if (selectedFile == null || !selectedFile.getPath().startsWith(ServerFile.SERVER_SESSION_ROOT_FOLDER)) {
			return null;
		}
		String filename = selectedFile.getPath().substring(ServerFile.SERVER_SESSION_ROOT_FOLDER.length()+1);
		@SuppressWarnings("unchecked")
		List<DbSession> sessions = (List<DbSession>)fileChooser.getClientProperty("sessions");
		return sessionManager.findSessionWithName(sessions, filename);
	}

	private void update() {
		try {
			RemoteSessionChooserFactory.updateRemoteSessions(sessionManager, fileChooser);

			StorageEntryMessageListener reply = sessionManager.getStorageUsage();
			
			long quota = reply.getQuota();
			long quotaWarning = reply.getQuotaWarning();
			long diskUsage = reply.getStorageUsage();
			
			quotaBar.setMaximum((int) (quota / 1024 / 1024));
			quotaBar.setValue((int) (diskUsage / 1024 / 1024));
			String humanReadableDiskUsage = Strings.toHumanReadable(diskUsage, true, true);
			String humanReadableQuota = Strings.toHumanReadable(quota, true, true);
			quotaBar.setString("Disk usage: " + humanReadableDiskUsage + "B / " +humanReadableQuota + "B");
			if (diskUsage < quotaWarning) {
				quotaBar.setForeground(LOW_DISK_USAGE_COLOR);
			} else {
				quotaBar.setForeground(HIGH_DISK_USAGE_COLOR);
			}
			
			String quotaWarningString = Strings.toHumanReadable(quotaWarning, true, true).trim() + "B";
			String quotaString = Strings.toHumanReadable(quota, true, true).trim() + "B";

			String lowDiskUsage = ""
					+ "Store up to " + quotaWarningString + " as long as you want.";
			
			String highDiskUsage = ""
					+ "Store up to " + quotaString + ", but please remove "
					+ "your data when you aren't anymore actively working on it.";
			
//			String highDiskUsage = ""
//			+ "Store up to " + quotaString + " when you are "
//			+ "actively working on it. We'll remind you by email if you "
//			+ "haven't used your data lately.";
			
			lowDiskUsageText.setText(lowDiskUsage);
			highDiskUsageText.setText(highDiskUsage);
			
			selectedSessionChanged();

		} catch (MalformedURLException | JMSException | InterruptedException e) {
			app.reportException(e);
		}		
	}

	@Override
	public void propertyChange(PropertyChangeEvent evt) {
		
		if (JFileChooser.SELECTED_FILE_CHANGED_PROPERTY
				.equals(evt.getPropertyName())) {

			selectedSessionChanged();
		}
	}

	private void selectedSessionChanged() {
		try {
			DbSession session = getSelectedSession();
			String uuid = null;
			String name = null;
			if (session != null) {
				uuid = session.getDataId();
				name = session.getName();
			}
			boolean isExampleSession = name != null && name.startsWith(DerbyMetadataServer.DEFAULT_EXAMPLE_SESSION_FOLDER);
			removeButton.setEnabled(name != null && !isExampleSession);
			updateSessionPreview(uuid);
		} catch (MalformedURLException e) {
			app.reportException(e);
		}
	}

	private void updateSessionPreview(final String uuid) {
		previewExecutor.execute(new Runnable() {
			
			@Override
			public void run() {
				try {
					updateGui(" ");
					
					if (uuid != null) {
						String preview = " ";

						StorageEntryMessageListener reply = sessionManager.getStorageUsage();
						for (StorageEntry session : reply.getEntries()) {						
							if (uuid.equals(session.getID())) {
								preview += Strings.toHumanReadable(session.getSize(), true, true) + "B, ";
								preview += session.getDate();
							}
						}

						updateGui(preview);
					}
					
				} catch (Exception e) {
					app.reportException(e);
				}
			}

			private void updateGui(final String preview) {
				SwingUtilities.invokeLater(new Runnable() {
					
					@Override
					public void run() {
						previewLabel.setText(preview);
					}
				});				
			}
		});
	}
}
