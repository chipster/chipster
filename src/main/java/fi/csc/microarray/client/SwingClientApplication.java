package fi.csc.microarray.client;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Component;
import java.awt.KeyEventDispatcher;
import java.awt.KeyboardFocusManager;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.lang.Thread.UncaughtExceptionHandler;
import java.net.MalformedURLException;
import java.net.URL;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.CertificateException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;
import javax.net.ssl.SSLHandshakeException;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRootPane;
import javax.swing.JSplitPane;
import javax.swing.KeyStroke;
import javax.swing.LookAndFeel;
import javax.swing.SwingUtilities;
import javax.swing.Timer;
import javax.swing.UIDefaults;
import javax.swing.UIManager;
import javax.swing.UnsupportedLookAndFeelException;
import javax.swing.filechooser.FileFilter;
import javax.swing.plaf.ColorUIResource;
import javax.swing.plaf.FontUIResource;

import org.apache.log4j.Logger;

import com.jgoodies.looks.HeaderStyle;
import com.jgoodies.looks.Options;
import com.jgoodies.looks.plastic.Plastic3DLookAndFeel;
import com.jgoodies.looks.plastic.PlasticTheme;
import com.jgoodies.looks.plastic.theme.ExperienceBlue;
import com.jgoodies.uif_lite.panel.SimpleInternalFrame;

import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.dataimport.ImportScreen;
import fi.csc.microarray.client.dataimport.ImportSession;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.dataview.GraphPanel;
import fi.csc.microarray.client.dataview.TreePanel;
import fi.csc.microarray.client.dialog.ChipsterDialog;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.PluginButton;
import fi.csc.microarray.client.dialog.DialogInfo;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.dialog.DialogInfo.Type;
import fi.csc.microarray.client.dialog.ErrorDialogUtils;
import fi.csc.microarray.client.dialog.ImportSettingsAccessory;
import fi.csc.microarray.client.dialog.SessionRestoreDialog;
import fi.csc.microarray.client.dialog.URLImportDialog;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.ToolPanel;
import fi.csc.microarray.client.screen.ChildScreenPool;
import fi.csc.microarray.client.screen.HistoryScreen;
import fi.csc.microarray.client.screen.Screen;
import fi.csc.microarray.client.screen.ShowSourceScreen;
import fi.csc.microarray.client.screen.TaskManagerScreen;
import fi.csc.microarray.client.serverfiles.ServerFile;
import fi.csc.microarray.client.session.RemoteSessionChooserFactory;
import fi.csc.microarray.client.session.UserSession;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.methods.DataDetails;
import fi.csc.microarray.client.waiting.WaitGlassPane;
import fi.csc.microarray.client.workflow.WorkflowManager;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.exception.ErrorReportAsException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.module.basic.BasicModule.VisualisationMethods;
import fi.csc.microarray.module.chipster.KielipankkiModule;
import fi.csc.microarray.util.BrowserLauncher;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.GeneralFileFilter;
import fi.csc.microarray.util.SplashScreen;
import fi.csc.microarray.util.Strings;

/**
 * This class adds all GUI and Swing specific content to client functionality.
 * 
 * @author Aleksi Kallio, Janne Käki
 * 
 */
public class SwingClientApplication extends ClientApplication {
	
	public static enum SessionSavingMethod {
		LEAVE_DATA_AS_IT_IS,
		INCLUDE_DATA_INTO_ZIP,
		UPLOAD_DATA_TO_SERVER;
	}

	private static final int METADATA_FETCH_TIMEOUT_SECONDS = 15;
	private static final long SLOW_VISUALISATION_LIMIT = 5 * 1000;
	private static final long VERY_SLOW_VISUALISATION_LIMIT = 20 * 1000;

	/**
	 * Logger for this class
	 */
	private static Logger logger;

	private JFrame mainFrame = null;
	private JPanel rightSideViewChanger = null;
	private JSplitPane rightSplit = null;
	private JSplitPane leftSplit = null;
	private JSplitPane mainSplit = null;
	private MicroarrayMenuBar menuBar;
	private TaskManagerScreen taskManagerScreen;
	private StatusBar statusBar;

	private SimpleInternalFrame treeFrame;
	private SimpleInternalFrame graphFrame;
	private SimpleInternalFrame operationsFrame;
	private JPanel visualisationArea;

	private ChildScreenPool childScreens;
	private TreePanel tree;
	private GraphPanel graphPanel;
	private ToolPanel toolPanel;
	private VisualisationFrameManager visualisationFrameManager;
	private HistoryScreen historyScreen;

	private SplashScreen splashScreen;
	private ClientListener clientListener;
	private WaitGlassPane waitPanel = new WaitGlassPane();
	
	private static float fontSize = VisualConstants.DEFAULT_FONT_SIZE;

	private JFileChooser importExportFileChooser;
	private JFileChooser workflowFileChooser;	
	
	/**
	 * Listens to jobExecutor's state in general.
	 */
	private PropertyChangeListener jobExecutorChangeListener = new PropertyChangeListener() {	
		public void propertyChange(PropertyChangeEvent evt) {
			taskCountChanged((Integer)evt.getNewValue(), true);
			logger.debug("JobExecutor property changed event: " + evt.getPropertyName() + ": " + (Integer)evt.getNewValue());
		}		
	};

	public SwingClientApplication(ClientListener clientListener, AuthenticationRequestListener overridingARL, String module)
	        throws MicroarrayException, IOException, IllegalConfigurationException {

		super(overridingARL);
		
		// this had to be delayed as logging is not available before loading configuration
		logger = Logger.getLogger(SwingClientApplication.class);
		
		Thread.currentThread().setUncaughtExceptionHandler(new UncaughtExceptionHandler() {				
			@Override
			public void uncaughtException(Thread t, Throwable e) {
				// we'll always output these to console and log for traceability and
				// easier IDE navigation
				e.printStackTrace();
				logger.error("Uncaught exception in thread " + t.getName(), e);
			}
		});
		
		if (!SwingUtilities.isEventDispatchThread()) {
			logger.error(new MicroarrayException("SwingClientApplication was created outside the Event Dispatch Thread."));
			System.exit(1);
		}
		
		this.clientListener = clientListener;

        // set the module that user wants to load
        this.requestedModule = module;

        // show splash screen
		if (KielipankkiModule.class.getName().equals(module)) {
			splashScreen = new SplashScreen(VisualConstants.getIcon(VisualConstants.SPLASH_SCREEN_KIELIPANKKI));
		} else {
			splashScreen = new SplashScreen(VisualConstants.getIcon(VisualConstants.SPLASH_SCREEN));
		}
		reportInitialisationThreadSafely("Initialising " + ApplicationConstants.TITLE, true);

		// try to initialise and handle exceptions gracefully
		
		/*Wait descriptions in another thread and let EDT continue.
		 * We can't wait in EDT, because otherwise there is a deadlock: LoginWindow
		 * waits for EDT to get free and EDT waits for LoginWindow to get done with authentication. 
		 */
		Thread t = new Thread(new Runnable() {
			
			@Override
			public void run() {
				try {
					initialiseApplication(false);

				} catch (Exception e) {
					if (Exceptions.isCausedBy(e, new SSLHandshakeException(""))) {
						reportTruststoreError(e);						
					} else {				
						reportInitalisationErrorThreadSafely(e);
					}
				}
			}

		});
		t.start();
	}

	private void reportTruststoreError(Exception e) {
		String truststore = null;
		try {
			truststore = JMSMessagingEndpoint.getClientTruststore();
		} catch (NoSuchAlgorithmException | CertificateException | KeyStoreException | IOException e1) {
			reportInitalisationErrorThreadSafely(new MicroarrayException("could not read trusstore configuration", e));
		}
			
		String title;
		String msg;
		PluginButton button;

		if (truststore == null) {
			title = "Server's identity cannot be verified";
			msg = "Server's administrator "
					+ "must install a proper certificate signed by CA or configure "
					+ "the server to use self-signed certificates.";
			button = null;
		} else {
			title = "Server's identity has changed";
			msg = "Click 'Close' to keep the old certificate. Server's identity will be verified "
					+ "with this certificate also next time."
					+ "\n\n"
					+ "Click 'Trust this server' to update the certificate. "
					+ "Update certificate only in a network that you trust. Chipster must be "
					+ "restarted manually after this.";
			
			final File truststoreFile = new File(truststore);
			
			button = new PluginButton() {

				@Override
				public void actionPerformed() {																		
					truststoreFile.delete();
					// now that the file is missing, getClientTruststore() will download it without complaints
					try {
						JMSMessagingEndpoint.getClientTruststore();
					} catch (NoSuchAlgorithmException | CertificateException | KeyStoreException | IOException e) {
						reportInitalisationErrorThreadSafely(new MicroarrayException("could not download trusstore", e));
					}
				}

				@Override
				public String getText() {
					return "Trust this server";
				}
			};
		}						
		showDialog(title, msg, Exceptions.getStackTrace(e), Severity.WARNING, true, DetailsVisibility.DETAILS_HIDDEN, button);
		quitImmediately();		
	}
	
	private void reportInitalisationErrorThreadSafely(final Exception e) {
		SwingUtilities.invokeLater(new Runnable() {
			
			@Override
			public void run() {

				showDialog("Starting Chipster failed.", "There could be a problem with the network connection, or the remote services could be down. " +
						"Please see the details below for more information about the problem.\n\n", 
						Exceptions.getStackTrace(e), Severity.ERROR, false);
				splashScreen.close();
				logger.error(e);			
			}
		});
	}

	public void reportInitialisationThreadSafely(final String report, final boolean newline) {
		SwingUtilities.invokeLater(new Runnable() {
			
			@Override
			public void run() {
				if (newline) {
					splashScreen.writeLine(report);
				} else {
					splashScreen.write(report);
				}			
			}
		});
	}

	public void initialiseGUIThreadSafely(final File mostRecentDeadTempDirectory) {
		SwingUtilities.invokeLater(new Runnable() {
			@Override
			public void run() {
				try {
					initialiseGUI(mostRecentDeadTempDirectory);
				} catch (MicroarrayException | IOException e) {
					reportInitalisationErrorThreadSafely(e);
				}
			}
		});
	}
	
	private void initialiseGUI(File mostRecentDeadTempDirectory) throws MicroarrayException, IOException {

		// assert state of initialisation
		try {
			definitionsInitialisedLatch.await(METADATA_FETCH_TIMEOUT_SECONDS, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
		if (toolModules == null) {
			throw new MicroarrayException("getting tool descriptions failed");
		}
		

		// initialize the main frame
		this.mainFrame = new JFrame();
		updateWindowTitle();
		childScreens = new ChildScreenPool(mainFrame);
		Frames frames = new Frames(mainFrame);
		Session.getSession().setFrames(frames);

		// set look'n'feel
		setPlastic3DLookAndFeel(mainFrame);

		// set location
		mainFrame.setLocationByPlatform(true);

		// initialise joblist popup menu
		this.taskManagerScreen = this.getTaskManagerScreen(); // call this method before getStatusBar to avoid null pointer exception

		// initialise child screens
		historyScreen = new HistoryScreen();

		childScreens.put("History", historyScreen);
		childScreens.put("ShowSource", new ShowSourceScreen());
		childScreens.put("Import", new ImportScreen());
		childScreens.put("TaskList", taskManagerScreen);

		// create operation panel using metadata
		try {
			toolPanel = new ToolPanel(toolModules);
		} catch (ParseException e) {
			logger.error("SADL parse failed", e);
			throw new MicroarrayException(e);
		}

		operationsFrame = getOperationsFrame();

		visualisationFrameManager = new VisualisationFrameManager();
		visualisationArea = getVisualisationFrameManager().getFramesPanel();

		rightSplit = new JSplitPane(JSplitPane.VERTICAL_SPLIT, operationsFrame, visualisationArea);
		rightSplit.setDividerLocation(VisualConstants.TREE_PANEL_HEIGHT);
		rightSplit.setResizeWeight(0.1);

		/* Initialize tree and graph */
		
		//moved to getTreeFrame
//		this.tree = new TreePanel(manager.getRootFolder());
		this.graphPanel = new GraphPanel();

		treeFrame = getTreeFrame();
		graphFrame = getGraphFrame();

		leftSplit = new JSplitPane(JSplitPane.VERTICAL_SPLIT, treeFrame, graphFrame);
		leftSplit.setDividerLocation(VisualConstants.TREE_PANEL_HEIGHT);
		leftSplit.setResizeWeight(0.1);


		rightSideViewChanger = new JPanel(new BorderLayout());

		if (this.isStandalone()) {
			rightSideViewChanger.add(visualisationArea, BorderLayout.CENTER);
		} else {	
			rightSideViewChanger.add(rightSplit, BorderLayout.CENTER);
		}
		rightSideViewChanger.setBorder(BorderFactory.createEmptyBorder());

		// construct the whole main content pane
		mainSplit = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftSplit, rightSideViewChanger);
		mainSplit.setDividerLocation(VisualConstants.LEFT_PANEL_WIDTH);
		mainSplit.setResizeWeight(0.1);

		// add menus
		menuBar = new MicroarrayMenuBar(this);
		
		// create status bar
		statusBar = new StatusBar(this);

		// sets 3D-look (JGoodies property)
		menuBar.putClientProperty(Options.HEADER_STYLE_KEY, HeaderStyle.SINGLE);

		// put everything together
		mainFrame.getContentPane().setLayout(new BorderLayout());
		mainFrame.getContentPane().add(mainSplit, BorderLayout.CENTER);
		mainFrame.getContentPane().add(statusBar.getStatusPanel(), BorderLayout.SOUTH);
		mainFrame.setJMenuBar(menuBar);
		menuBar.updateMenuStatus();

		// add glass wait panel
		JRootPane rootPane = SwingUtilities.getRootPane(mainFrame);
		rootPane.setGlassPane(waitPanel);

		// add shutdown listener
		mainFrame.addWindowListener(new java.awt.event.WindowAdapter() {
			public void windowClosing(java.awt.event.WindowEvent e) {
				quit();
			}
		});

		// make window visible
		mainFrame.setIconImage(VisualConstants.getIcon(VisualConstants.APPLICATION_ICON).getImage());
		mainFrame.pack();
		mainFrame.setExtendedState(JFrame.MAXIMIZED_BOTH);
		mainFrame.setVisible(true);
		mainFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		
		// start listening to job events
		taskExecutor.addChangeListener(jobExecutorChangeListener);
		
		// Start checking amount of free memory 
		final Timer memoryCheckTimer = new Timer(MEMORY_CHECK_INTERVAL, new ActionListener() {
			public void actionPerformed(ActionEvent e) {
				checkFreeMemory();
			}
		});
		memoryCheckTimer.setCoalesce(true);
		memoryCheckTimer.setRepeats(true);
		memoryCheckTimer.setInitialDelay(0);
		memoryCheckTimer.start();
		
		// it's alive!
		super.setEventsEnabled(true);
		manager.setEventsEnabled(true);

		// hide splashscreen
		splashScreen.close();

		// notify listener
		if (clientListener != null) {
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					clientListener.onSuccessfulInitialisation();
				}
			});
		}

		// final touches
		updateFocusTraversal();
		restoreDefaultView();
		enableKeyboardShortcuts();		
		setVisualisationMethodToDefault();
		
		// check for warnings generated at earlier non-GUI stages
		if (!getInitialisationWarnings().isEmpty()) {
			showDialog("Missing tools", "Some tools were not loaded during the startup and they will not be shown in the tool lists. This suggests that there are technical problems with the server you are connected to.", getInitialisationWarnings(), Severity.WARNING, true);
		}
		
		// check for session restore need
		if (mostRecentDeadTempDirectory != null) {
			
			File sessionFile = UserSession.findBackupFile(mostRecentDeadTempDirectory, false);
			new SessionRestoreDialog(this, sessionFile).setVisible(true);
		}
	}
	
	/**
	 * Set default visualisation method (DataDetails)
	 */
	public void setVisualisationMethodToDefault() {
		setVisualisationMethod(null, null, getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
	}

	protected void showDebugDialog(int type) {
		
		switch (type) {
		default:
			String dump = getSessionManager().printSession();
			ChipsterDialog.showDialog(this, new DialogInfo(Severity.INFO, "Session URL dump", "See details for URL's of data beans.", dump), DetailsVisibility.DETAILS_ALWAYS_VISIBLE, false);
		}
	}

	public void updateFocusTraversal() {
		Vector<Component> order = new Vector<Component>();
		order.addAll(visualisationFrameManager.getFocusComponents());
		order.addAll(tree.getFocusComponents());
		order.addAll(toolPanel.getFocusComponents());

		getMainFrame().setFocusTraversalPolicy(new ClientFocusTraversalPolicy(order));

	}

	
	private String windowTitleJobPrefix = null;
	private String windowTitleBlockingPrefix = null;
	private JFileChooser localSessionFileChooser;

	public void updateWindowTitleJobCount(Integer jobCount) {
		windowTitleJobPrefix = jobCount > 0 ? jobCount + " tasks / " : null;
		updateWindowTitle();
	}
	
	public void updateWindowTitleBlockingState(String operation) {
		if (operation != null) {
			windowTitleBlockingPrefix = Strings.startWithUppercase(operation) + " / ";
			
		} else {
			windowTitleBlockingPrefix = null;
		}	
		updateWindowTitle();
	}
	
	public void updateWindowTitle() {
		if (windowTitleBlockingPrefix != null) {
			this.mainFrame.setTitle(windowTitleBlockingPrefix + Session.getSession().getPrimaryModule().getDisplayName() + " " + ApplicationConstants.VERSION);
			
		} else if (windowTitleJobPrefix != null) {
			this.mainFrame.setTitle(windowTitleJobPrefix + Session.getSession().getPrimaryModule().getDisplayName() + " " +  ApplicationConstants.VERSION);
			
		} else {
			this.mainFrame.setTitle(Session.getSession().getPrimaryModule().getDisplayName() + " " + ApplicationConstants.VERSION);
		}

	}

	public SimpleInternalFrame getOperationsFrame() {
		if (operationsFrame == null) {
			operationsFrame = new SimpleInternalFrame("Analysis tools");
			operationsFrame.setContent(toolPanel);
		}
		return operationsFrame;
	}

	public SimpleInternalFrame getGraphFrame() throws MicroarrayException {
		if (graphFrame == null) {
			graphFrame = new SimpleInternalFrame("Workflow", graphPanel.getButtonToolBar(), graphPanel.getScroller());
		}
		return graphFrame;
	}

	public SimpleInternalFrame getTreeFrame() throws MicroarrayException {
		if (treeFrame == null) {
			treeFrame = new SimpleInternalFrame("Datasets");
						
			CardLayout cardLayout = new CardLayout();
			JPanel cardPanel = new JPanel(cardLayout);
			
			QuickLinkPanel linkPanel = new QuickLinkPanel();
			this.tree = new TreePanel(manager.getRootFolder(), cardPanel, cardLayout);
						
			cardPanel.add(linkPanel, "LINKS");
			cardPanel.add(tree, "TREE");
			
			treeFrame.add(cardPanel);
			cardLayout.first(cardPanel);
			
			return treeFrame;

		} else {
			return treeFrame;
		}
	}



	/**
	 * ExperienceBlue is very nice color theme, but it has ugly orange borders
	 * around focused components. This helper class customizes the theme a bit.
	 * 
	 * @author mkoski
	 */
	private static class CustomExperienceBlue extends ExperienceBlue {

		/**
		 * Removes the ugly orange focus color
		 */
		@Override
		public ColorUIResource getFocusColor() {
			return new ColorUIResource(VisualConstants.PLASTIC3D_FOCUS_COLOR);
		}

		@Override
		public FontUIResource getControlTextFont() {
			return new FontUIResource(super.getControlTextFont().deriveFont(fontSize));
		}

		@Override
		public FontUIResource getTitleTextFont() {
			return new FontUIResource(super.getTitleTextFont().deriveFont(fontSize));
		}

		@Override
		public FontUIResource getSubTextFont() {
			return new FontUIResource(super.getSubTextFont().deriveFont(fontSize));
		}

		@Override
		public FontUIResource getSystemTextFont() {
			return new FontUIResource(super.getSystemTextFont().deriveFont(fontSize));
		}

		@Override
		public FontUIResource getUserTextFont() {
			return new FontUIResource(super.getUserTextFont().deriveFont(fontSize));
		}

		@Override
		public FontUIResource getWindowTitleFont() {
			return new FontUIResource(super.getWindowTitleFont().deriveFont(fontSize));
		}

		@Override
		public FontUIResource getMenuTextFont() {
			return new FontUIResource(super.getMenuTextFont().deriveFont(fontSize));
		}
	}

	private static PlasticTheme theme;
	private static LookAndFeel lnf;

	/**
	 * This method sets applications appearance to new level by using Plastic3D
	 * Look And Feel and ExperienceBlue color theme.
	 * 
	 * The method sets look and feel and updates UIDefault values. After that it
	 * updates all the components to the new appearance.
	 * 
	 * @param componentTreeRoot
	 *            root component of component tree (frame component in most
	 *            cases)
	 * 
	 */
	public static void setPlastic3DLookAndFeel(Component componentTreeRoot) {

		// Look and Feel change must be done only once, otherwise the custom
		// values
		// from the components get overridden. However this method is called
		// several
		// times for different child windows, but these ifs do the job.
		if (theme == null) {
			theme = new CustomExperienceBlue();
			Plastic3DLookAndFeel.setPlasticTheme(theme);
		}
		if (lnf == null) {
			lnf = new Plastic3DLookAndFeel();

			try {
				UIManager.installLookAndFeel("Plastic3D", lnf.getClass().getName());
				UIManager.setLookAndFeel(lnf);
			} catch (UnsupportedLookAndFeelException e) {
				e.printStackTrace();
				logger.error(Exceptions.getStackTrace(e));
			}
		}

		// set the UI defaults
		UIDefaults defaults = UIManager.getDefaults();

		// Done several times, but luckily map rejects duplicates
		defaults.putAll(VisualConstants.getUIDefaults());

		SwingUtilities.updateComponentTreeUI(componentTreeRoot);
	}

	public void setFontSize(float size) {
		fontSize = size;

		SwingUtilities.updateComponentTreeUI(mainFrame);

		fixFileChooserFontSize(importExportFileChooser);

		// Running progressBar don't care about updateUI
		statusBar.setFontSize(fontSize);

		Iterator<Screen> iter = childScreens.getScreenIterator();

		while (iter.hasNext()) {
			Screen screen = iter.next();
			SwingUtilities.updateComponentTreeUI(screen.getFrame());
		}
	}

	public float getFontSize() {
		return fontSize;
	}

	public JFrame getMainFrame() {
		return mainFrame;
	}




	/**
	 * Sets the folder with given name selected or creates a new folder if it
	 * doesn't exist yet. Selected folder is used as a place for imported files.
	 * If folder name is null, root folder is returned
	 * 
	 * @param folderName
	 *            folder name
	 * @return just created or previously existing folder. If given folder name
	 *         is <code>null</code>, root folder is returned
	 */
	public DataFolder initializeFolderForImport(String folderName) {

		DataFolder root = manager.getRootFolder();

		if (folderName == null || folderName.isEmpty()) {
			logger.debug("initializing for import " + folderName + ": is null => using root");
			return root;

		} else if (ImportUtils.getFolderNames(false).contains(folderName)) {

			logger.debug("initializing for import " + folderName + ": exists already");
			DataFolder folderToSelect;
			if (folderName.equals(root.getName())) {
				folderToSelect = root;
			} else {
				folderToSelect = root.getChildFolder(folderName);
			}

			getSelectionManager().selectSingle(folderToSelect, this);
			return folderToSelect;

		} else {

			logger.debug("initializing for import " + folderName + ": creating new ");
			DataFolder folder = getDataManager().createFolder(root, folderName);
			getSelectionManager().selectSingle(folder, this);
			return folder;
		}
	}

	public void importGroup(final Collection<ImportItem> datas, final String folderName) {

		runBlockingTask("importing files", new Runnable() {

			public void run() {
				importGroupAndWait(datas, folderName);
			}
		});
	}

	/**
	 * Convenience method for showing exception dialogs with hand written
	 * messages and stack traces hidden first.
	 * 
	 * @see #showDialog(String, String, String, Severity)
	 */

	public void showErrorDialog(String title, Exception error) {
		showDialog(title, null, Exceptions.getStackTrace(error), Severity.ERROR, false);
	}

	/**
	 * @see #showDialog(String, String, String, Severity)
	 */
	public void showDialog(String title, Severity severity, boolean modal) {
		showDialog(title, null, null, severity, modal);
	}

	public void showDialog(String title, String message, String details, Severity severity, boolean modal) {
		showDialog(title, message, details, severity, modal, ChipsterDialog.DetailsVisibility.DETAILS_HIDDEN, null);
	}

	/**
	 * Shows a modal Chipster-styled dialog box to user. Use only when you need
	 * user's immediate attention.
	 * 
	 * @param message
	 *            message that is always shown to user
	 * @param details
	 *            information that is shown in text box that if first closed
	 * @param severity
	 *            severity level, affects icon choice
	 */
	public void showDialog(String title, String message, String details, Severity severity, boolean modal, ChipsterDialog.DetailsVisibility detailsVisibility, PluginButton button) {
		DialogInfo dialogInfo = new DialogInfo(severity, title, message, details);
		ChipsterDialog.showDialog(this, dialogInfo, detailsVisibility, modal, null, button);
	}

	public void showDialog(String title, String message, String details, Severity severity, boolean modal, ChipsterDialog.DetailsVisibility detailsVisibility, PluginButton button, boolean feedbackEnabled) {
		DialogInfo dialogInfo = new DialogInfo(severity, title, message, details);
		dialogInfo.setFeedbackVisible(feedbackEnabled);
		ChipsterDialog.showDialog(this, dialogInfo, detailsVisibility, modal, null, button);
	}

	public File saveWorkflow() {

		try {
			JFileChooser fileChooser = this.getWorkflowFileChooser();
			int ret = fileChooser.showSaveDialog(this.getMainFrame());
			if (ret == JFileChooser.APPROVE_OPTION) {
				File selected = fileChooser.getSelectedFile();
				File newFile = selected.getName().endsWith(WorkflowManager.SCRIPT_EXTENSION) ? selected : new File(selected.getCanonicalPath() + "." + WorkflowManager.SCRIPT_EXTENSION);

				super.saveWorkflow(newFile);
				menuBar.addRecentWorkflow(newFile.getName(), Files.toUrl(newFile));
				menuBar.updateMenuStatus();
				return newFile;
			}
			menuBar.updateMenuStatus();

		} catch (IOException e) {
			reportException(e);
		}
		return null;
	}

	public void runWorkflow(final URL workflowScript, final boolean runForEach) {
		Thread thread = new Thread(new Runnable() {

			@Override
			public void run() {
				
				SwingClientApplication.super.runWorkflow(workflowScript, runForEach);
			}
		});
		thread.start();
	}

	public File openWorkflow(boolean runForEach) {

		try {
			JFileChooser fileChooser = this.getWorkflowFileChooser();
			int ret = fileChooser.showOpenDialog(this.getMainFrame());
			if (ret == JFileChooser.APPROVE_OPTION) {
				runWorkflow(fileChooser.getSelectedFile().toURI().toURL(), runForEach);

				menuBar.updateMenuStatus();
				return fileChooser.getSelectedFile();
			} else {
				menuBar.updateMenuStatus();
				return null;
			}
			
		} catch (MalformedURLException e) {
			reportException(e);
			return null;
		}
	}

	public void showHistoryScreenFor(DataBean data) {
		historyScreen.setData(data);
		childScreens.show("History", true);
	}

	public void restoreDefaultView() {
		leftSplit.setDividerLocation(VisualConstants.TREE_PANEL_HEIGHT);
		mainSplit.setDividerLocation(VisualConstants.LEFT_PANEL_WIDTH);
		rightSplit.setDividerLocation(VisualConstants.TREE_PANEL_HEIGHT);
		mainSplit.validate();
	}

	public void reportTaskError(Task task) throws MicroarrayException {
		String title;
		String message;
		
		boolean userFixable = task.getState() == State.FAILED_USER_ERROR && task.getErrorMessage() != null && !task.getErrorMessage().equals("");
		
		// user-friendly message
		if (userFixable) {
			title = task.getErrorMessage();
			message = task.getFullName() + " was stopped. ";
		} 
		
		// generic message
		else {
			title = task.getFullName() + " did not complete successfully. ";
			message = "You may have used a tool or parameters which are unsuitable for the selected dataset, or " + "there might be a bug in the analysis tool itself.\n\n" + "The details below may provide hints about the problem.";
		}		

		// details
		String details = "";
		if (task.getErrorMessage() != null) {
			details += task.getErrorMessage() + "\n\n";
		}
		
		if (task.getStateDetail() != null) {
			details += task.getStateDetail() + "\n\n";
		}
		
		details += "----------------------------------------------------------------------\n";
		
		if (task.getScreenOutput() != null) {
			details += task.getScreenOutput();
		}
		
		DetailsVisibility detailsVisibility = userFixable ? DetailsVisibility.DETAILS_ALWAYS_HIDDEN : DetailsVisibility.DETAILS_HIDDEN;
		
		// show dialog
		DialogInfo dialogInfo = new DialogInfo(Severity.INFO, title, message, details);
		if (!userFixable) {
			dialogInfo.setFeedbackVisible(true);
		}
		
		ChipsterDialog.showDialog(this, dialogInfo, detailsVisibility, false);
	}
	
	public void reportExceptionThreadSafely(final Exception e) {
		SwingUtilities.invokeLater(new Runnable() {			
			@Override
			public void run() {
				reportException(e);
			}
		});
	}

	public void reportException(Exception e) {

		// collect error information to dialogInfo
		DialogInfo dialogInfo = new DialogInfo(Severity.ERROR,
		        "An error has occurred and the action was not performed successfully.",
		        "If problem persists, please check that your data is valid. For more information open the details panel below.", null);
		dialogInfo.setFeedbackVisible(true);

		// exception has extra info
		if (e instanceof ErrorReportAsException) {

			// exception has everything we need
			ErrorReportAsException report = (ErrorReportAsException) e;
			dialogInfo.setTitle(report.getTitle());
			dialogInfo.setMessage(report.getMessage());
			dialogInfo.setDetails(report.getDetails());

		} else if (e instanceof MicroarrayException) {

			// exception can be scavenged for details
			MicroarrayException me = (MicroarrayException) e;
			String details = "";
			
			if (e.getMessage() != null) {
				details += e.getMessage() + "\n";
			}

			if (ErrorDialogUtils.getMessage(me.getExtraInfo()) != null) {
				details += ErrorDialogUtils.getMessage(me.getExtraInfo()) + "\n";
			}

			if (ErrorDialogUtils.getScreenOutput(me.getExtraInfo()) != null) {
				details += ErrorDialogUtils.getScreenOutput(me.getExtraInfo());
			}

			dialogInfo.setDetails(details);

		} else {
			String details = "";
			
			// use stack trace as details
			details += Exceptions.getStackTrace(e);
			
			if (e.getCause() != null) {
				details += "Caused by: ";
				details += Exceptions.getStackTrace(e.getCause());				
			}
			dialogInfo.setDetails(details);
		}

		// show dialog
		ChipsterDialog.showDialog(this, dialogInfo, ChipsterDialog.DetailsVisibility.DETAILS_HIDDEN, false);

		// we'll always output these to console and log for traceability and
		// easier IDE navigation
		e.printStackTrace();
		if (logger != null) {
			logger.error(e.getMessage(), e);
			e.printStackTrace();
		}
	}

	public void onException(JMSException e) {
		reportException(e);
	}

	public void showChildScreen(String name, boolean packed) {
		childScreens.show(name, packed);
	}

	public void showMaximisedVisualisation(boolean maximised) {
		if (maximised) {
			mainSplit.setDividerLocation(0);
			rightSplit.setDividerLocation(0);
		} else {
			restoreDefaultView();
		}
		rightSideViewChanger.validate();
	}

	public void showPopupMenuFor(MouseEvent e, DataItem data) {
		List<DataItem> datas = new ArrayList<DataItem>();
		datas.add(data);
		showPopupMenuFor(e, datas);
	}

	public void showPopupMenuFor(MouseEvent e, List<DataItem> datas) {
		ClientContextMenu popup = new ClientContextMenu(this);
		popup.setOptionsFor(datas);
		popup.show(e.getComponent(), e.getX(), e.getY());
	}

	/**
	 * Is the selected databeans possible to visualise
	 */
	public boolean isSelectedDataVisualisable() {
		return getDefaultVisualisationForSelection() != null;
	}



	public void visualiseWithBestMethod(FrameType target) {
		setVisualisationMethod(getDefaultVisualisationForSelection(), null, getSelectionManager().getSelectedDataBeans(), target);
	}

	public void deleteDatas(DataItem... datas) {

		// check that we have something to delete
		if (datas.length == 0) {
			return; // no selection, do nothing
		}

		// choose confirm dialog message
		JLabel confirmMessage = null;
		if (datas.length > 1) {
			confirmMessage = new JLabel("Really delete " + datas.length + " data items?");
			
		} else if (datas[0] instanceof DataFolder) {
			confirmMessage = new JLabel("Really delete " + datas[0].getName() + " and all of its contents?");
			
		} else if (datas[0] instanceof DataBean) {
			confirmMessage = new JLabel("Really delete " + datas[0].getName() + " ?");
			
		} else {
			throw new IllegalArgumentException("datas is illegal");
		}
		
		// confirm delete
		if (JOptionPane.showConfirmDialog(this.mainFrame, confirmMessage, "Delete items", JOptionPane.YES_NO_OPTION) != JOptionPane.YES_OPTION) {
			return; // deletion was not confirmed
		}

		// delete actually
		super.deleteDatasWithoutConfirming(datas);
	}


	public void fixFileChooserFontSize(JFileChooser fileChooser) {
		// Some special care has to be taken to get the fileChooser list update
		// also
		if (fileChooser != null) {
			UIManager.put("FileChooser.listFont", UIManager.getFont("List.font").deriveFont(fontSize));
			SwingUtilities.updateComponentTreeUI(fileChooser);
		}
	}

	public void openFileImport() throws MicroarrayException, IOException {
		JFileChooser fc = getImportExportFileChooser();
		fc.setMultiSelectionEnabled(true);

		ImportSettingsAccessory access = (ImportSettingsAccessory) importExportFileChooser.getAccessory();
		access.setDefaults();
		int ret = fc.showOpenDialog(getMainFrame());
		if (ret == JFileChooser.APPROVE_OPTION) {
			List<Object> files = new LinkedList<Object>();

			for (File file : fc.getSelectedFiles()) {
				files.add(file);
			}

			ImportSession importSession = new ImportSession(ImportSession.Source.FILE, files, access.getImportFolder(), access.skipActionChooser());
			ImportUtils.executeImport(importSession);
		}
	}

	public void openDirectoryImportDialog() {
		JFileChooser fc = getImportExportDirChooser();
		fc.setSelectedFile(new File(""));
		int ret = fc.showOpenDialog(this.getMainFrame());
		if (ret == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fc.getSelectedFile();
			super.importWholeDirectory(selectedFile);
		}
	}

	public void openURLImport() throws MicroarrayException, IOException {
		URLImportDialog urlImportDlg = new URLImportDialog(this);
		URL selectedURL = urlImportDlg.getSelectedURL();
		String importFolder = urlImportDlg.getSelectedFolderName();
		if (selectedURL != null) {
			ImportUtils.executeImport(new ImportSession(ImportSession.Source.URL, new URL[] { selectedURL }, importFolder, urlImportDlg.isSkipSelected()));
		}
	}
	
	public void quit() {
				
		if (!killUploadingTasks()) {
			// user wanted to continue uploading, so we can't quit
			return;
		}

		// Check for unsaved changes
		int returnValue = JOptionPane.DEFAULT_OPTION;

		if (getSessionManager().hasUnsavedChanges()) {

			Object[] options = { "Save", "Close without saving", "Cancel" };

			returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), "Do you want the session to be saved to server before closing Chipster?", "Confirm close", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

			if (returnValue == 0) {
				try {
					if (getSessionManager().areCloudSessionsEnabled()) {
						saveSession(SessionSavingMethod.UPLOAD_DATA_TO_SERVER);
					} else {
						saveSession(SessionSavingMethod.INCLUDE_DATA_INTO_ZIP);
					}
					return;
				} catch (Exception exp) {
					this.showErrorDialog("Session saving failed", exp);
					return;
				}

				// cancel also if the dialog is closed from the corner x
			} else if (returnValue != 1) {
				return;
			}
		}
		
		quitImmediately();
	}
	
	/**
	 * Check if there are uploading jobs and ask if user wants to kill
	 * them. Returns true if the jobs were killed or there weren't any. Returns
	 * false if the killing was cancelled.
	 * 
	 * @return false if the action was cancelled
	 */
	private boolean killUploadingTasks() {
		// Check the running tasks
		int returnValue = JOptionPane.DEFAULT_OPTION;
		
		if (taskExecutor.getUploadingTaskCount() > 0) {
			String message = "";
			if (taskExecutor.getUploadingTaskCount() == 1) {
				message += "There is a task uploading input files.  Are you sure you want to cancel the task?";
			} else {
				message += "There are " + taskExecutor.getUploadingTaskCount() + " tasks uploading input files. " + "Are you sure you want to cancel these tasks?";
			}

			Object[] options = { "Cancel uploading tasks", "Continue uploading" };

			returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), message, "Confirm close", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

			if (returnValue == JOptionPane.YES_OPTION) {
				taskExecutor.killUploadingTasks();
			} else {
				return false;
			}
		}
		return true;
	}

	public void quitImmediately() {

		// hide immediately to look more reactive...
		if (childScreens != null) {
			childScreens.disposeAll();
		}
		if (mainFrame != null) {
			mainFrame.setVisible(false);
		}

		super.quit();

		// this closes the application
		if (mainFrame != null) {
			mainFrame.dispose();
		}
		System.exit(0);
	}

	private static ClientListener getShutdownListener() {
		ClientListener shutdownListener = new ClientListener() {
			public void onSuccessfulInitialisation() {
				// do nothing
			}
			public void onFailedInitialisation() {
				System.exit(1);
			}
		};
		return shutdownListener;
	}
	

	/**
	 * Starts Chipster client. Configuration (logging) should be initialised
	 * before calling this method.
	 */
	public static void start(String configURL, String module) throws IOException {

		try {
			DirectoryLayout.initialiseClientLayout(configURL);			

		} catch (IllegalConfigurationException e) {
			reportIllegalConfigurationException(e);
		}

		ClientListener shutdownListener = getShutdownListener();
		
		try {						
			new SwingClientApplication(shutdownListener, null, module);
			
		} catch (Throwable t) {
			t.printStackTrace();
			if (logger != null) {
				logger.error(Exceptions.getStackTrace(t));
			}
		}

	}
	
	public static void main(String[] args) throws IOException {
		start(null, "fi.csc.microarray.module.chipster.MicroarrayModule");
	}

	public static void reportIllegalConfigurationException(IllegalConfigurationException e) {
		DialogInfo dialogInfo = new DialogInfo(Severity.ERROR, "Illegal configuration", "Chipster could not start because the provided configuration file is illegal. Please contact your system administrator.", "Reason: " + e.getMessage());
		ChipsterDialog.showDialog(null, dialogInfo, DetailsVisibility.DETAILS_HIDDEN, true);
		throw new RuntimeException("configuration not compatible, will not start");
	}

	public void showSourceFor(String operationID) throws TaskException {
		childScreens.show("ShowSource", true, operationID);
	}

	/**
	 * Opens JFileChooser to ask location for the export, then calls
	 * exportDataDirectory to do the job.
	 * 
	 * @param data
	 * @throws MicroarrayException
	 * @throws IOException
	 */
	public void exportFolder(DataFolder data) throws MicroarrayException, IOException {

		File file = new File(data.getName().replace(" ", "_"));
		// file.mkdirs();

		JFileChooser fc = getImportExportDirChooser();
		fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		fc.setSelectedFile(file);

		logger.debug("Exporting File: " + fc.getSelectedFile().getAbsolutePath());

		int ret = fc.showSaveDialog(getMainFrame());

		if (ret == JFileChooser.APPROVE_OPTION) {
			File selected = new File(fc.getSelectedFile().getAbsolutePath() + File.separator + file.getName());
			selected.mkdirs();
			exportDataFolder(data, selected);
		}
	}

	/**
	 * Method to export recursively all descendants of the source
	 * 
	 */
	private void exportDataFolder(DataFolder source, File dir) throws IOException, MicroarrayException {
		for (DataItem child : source.getChildren()) {
			if (child instanceof DataBean) {
				String filename = createFilename((DataBean) child);
				filename = dir.getAbsolutePath() + File.separator + filename;
				logger.debug("Exporting dataBean " + child.getName() + " into " + filename);
				exportToFile((DataBean) child, new File(filename));

			} else if (child instanceof DataFolder) {
				logger.debug("Exporting dataFolder " + child.getName());
				String foldername = dir.getAbsolutePath() + File.separator + child.getName().replace(" ", "_");

				File file = new File(foldername);
				file.mkdir();
				exportDataFolder((DataFolder) child, file);
			}
		}
	}
	
	/**
	 * FIXME Better handling for existing file
	 * 
	 * @param data
	 * @param selectedFile
	 */
	protected void exportToFile(final DataBean data, final File selectedFile) {
		runBlockingTask("exporting file", new Runnable() {

			public void run() {
				exportToFileAndWait(data, selectedFile);
			}			
		});		
	}

	public void exportDataset(DataBean data) throws MicroarrayException, IOException {
		JFileChooser fc = this.getImportExportFileChooser();
		fc.setDialogType(JFileChooser.SAVE_DIALOG);
		String proposedName = createFilename(data);
		fc.setSelectedFile(new File(proposedName));
		int ret = fc.showSaveDialog(getMainFrame());
		if (ret == JFileChooser.APPROVE_OPTION) {
			File selectedFile = fc.getSelectedFile();
			exportToFile(data, selectedFile);
		}
	}

	private String createFilename(DataBean data) {
		String proposedName = data.getName().replace(" ", "_");
		/*
		 * TODO does not work with mime content types if
		 * (!proposedName.endsWith(data.getContentType())) { proposedName =
		 * proposedName + "." + data.getContentType(); }
		 */
		return proposedName;
	}

	@Override
	public void setVisualisationMethod(VisualisationMethod method, List<Variable> variables, List<DataBean> datas, FrameType target) {
		
		if (method == null) {
			boolean datasetsSelected = (datas != null && !datas.isEmpty());
			boolean datasetsExist = !getDataManager().databeans().isEmpty();
			
			if (datasetsSelected) {
				method = VisualisationMethods.DATA_DETAILS;
			} else if (datasetsExist) {
				method = VisualisationMethods.SESSION_DETAILS;
			} else {				
				method = VisualisationMethods.EMPTY;
			}				
			
			super.setVisualisationMethod(method, variables, datas, target);
			return;
		}

		long estimate = method.estimateDuration(datas);

		if (estimate > SLOW_VISUALISATION_LIMIT) {
			int returnValue = JOptionPane.DEFAULT_OPTION;

			String message = "";
			int severity;
			// Check the running tasks
			if (estimate > VERY_SLOW_VISUALISATION_LIMIT) {

				message += "Visualising the selected large dataset with this method might stop Chipster from responding. \n" + "If you choose to continue, it's recommended to save the session before visualising.";
				severity = JOptionPane.WARNING_MESSAGE;
			} else {
				message += "Are you sure you want to visualise large dataset, which may " + "take several seconds?";
				severity = JOptionPane.QUESTION_MESSAGE;
			}
			Object[] options = { "Cancel", "Visualise" };

			returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), message, "Cancel visualisation", JOptionPane.YES_NO_OPTION, severity, null, options, options[0]);

			if (returnValue == 1) {
				super.setVisualisationMethod(method, variables, datas, target);
			} else {
				return;
			}
		} else {
			super.setVisualisationMethod(method, variables, datas, target);
		}
	}

	public void exportSelectedItems() {
		try {

			for (DataItem item : getSelectionManager().getSelectedDataItems()) {
				if (item instanceof DataBean) {
					exportDataset((DataBean) item);
				} else if (item instanceof DataFolder) {
					exportFolder((DataFolder) item);
				}
			}
		} catch (Exception e) {
			reportException(e);
		}
	}

	private JFileChooser getImportExportDirChooser() {
		if (importExportFileChooser == null) {
			importExportFileChooser = ImportUtils.getFixedFileChooser();
		}

		importExportFileChooser.setAccessory(null);
		importExportFileChooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
		String description = "Data Directory";
		String[] extensions = {}; // only directories

		for (FileFilter filter : importExportFileChooser.getChoosableFileFilters()) {
			importExportFileChooser.removeChoosableFileFilter(filter);
		}

		importExportFileChooser.addChoosableFileFilter(new GeneralFileFilter(description, extensions));
		importExportFileChooser.setAcceptAllFileFilterUsed(false);

		fixFileChooserFontSize(importExportFileChooser);

		return importExportFileChooser;
	}

	/**
	 * This is public to give access to file filter from ImportUtils
	 * 
	 * @return The default JFileChooser with an assigned MicroarrayFileFilter.
	 */
	public JFileChooser getImportExportFileChooser() {
		if (importExportFileChooser == null) {
			importExportFileChooser = ImportUtils.getFixedFileChooser();

		}

		// remove previous filters, otherwise they duplicate
		for (FileFilter filter : importExportFileChooser.getChoosableFileFilters()) {
			importExportFileChooser.removeChoosableFileFilter(filter);
		}

		importExportFileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		importExportFileChooser.setAcceptAllFileFilterUsed(true);

		ImportSettingsAccessory access = new ImportSettingsAccessory(importExportFileChooser);
		importExportFileChooser.setAccessory(access);

		fixFileChooserFontSize(importExportFileChooser);

		return importExportFileChooser;
	}

	private JFileChooser getSessionFileChooser(boolean remote, boolean openExampleDir) throws MalformedURLException, JMSException, Exception {
		
		if (openExampleDir) {
			
			return new RemoteSessionChooserFactory(this).getExampleSessionChooser();
			
		} else if (remote) {
			
			return new RemoteSessionChooserFactory(this).getRemoteSessionChooser();

		} else {
			// don't create new local file choosers, because that would reset
			// the selected directory
			if (localSessionFileChooser == null) {

				localSessionFileChooser = ImportUtils.getFixedFileChooser();
				localSessionFileChooser.setFileFilter(new GeneralFileFilter("Chipster Session (." + UserSession.SESSION_FILE_EXTENSION + ")", new String[] { UserSession.SESSION_FILE_EXTENSION }));
				localSessionFileChooser.setAcceptAllFileFilterUsed(false);
				localSessionFileChooser.setMultiSelectionEnabled(false);
				fixFileChooserFontSize(localSessionFileChooser);
				localSessionFileChooser.setSelectedFile(new File("session." + UserSession.SESSION_FILE_EXTENSION));
			}
			return localSessionFileChooser;
		}

	}

	private JFileChooser getWorkflowFileChooser() {
		if (workflowFileChooser == null) {
			workflowFileChooser = ImportUtils.getFixedFileChooser(workflowManager.getScriptDirectory());

			workflowFileChooser.setFileFilter(WorkflowManager.FILE_FILTER);
			workflowFileChooser.setSelectedFile(new File("workflow.bsh"));
			workflowFileChooser.setAcceptAllFileFilterUsed(false);
			workflowFileChooser.setMultiSelectionEnabled(false);
		}

		fixFileChooserFontSize(workflowFileChooser);

		return workflowFileChooser;
	}

	/**
	 * Opens import tool directly
	 * 
	 * @param useSameSettings
	 */
	public void openImportTool(ImportSession importSession) {
		
		if (!importSession.isLocal()) {
			throw new IllegalArgumentException("only local files supported in import tool");
		}
		
		// fire up import tool
		ImportScreen importScreen = (ImportScreen) childScreens.get("Import");
		importScreen.setImportSession(importSession);
		importScreen.updateTable(false);
		childScreens.show("Import", false);
	}

	/**
	 * <p>Run a task in background thread and block GUI in a friendly way while the task is being
	 * run. <strong>Note!</strong> If the task modifies GUI, it must use SwingUtilities.invokeAndWait so that
	 * modifications are done in event dispatch thread.</p>
	 * 
	 * @param taskName
	 *            name of the task we will be running, like "importing" or
	 *            "loading"
	 * @param runnable
	 *            the task
	 */
	public void runBlockingTask(String taskName, final Runnable runnable) {
		
		Thread backgroundThread = new Thread(new Runnable() {
			public void run() {
				try {
					runnable.run();
				} catch (final Exception e) {
					SwingUtilities.invokeLater(new Runnable() {
						public void run() {
							reportException(e);
						}
					});
				} finally {
					waitPanel.stopWaiting();
					updateWindowTitleBlockingState(null);
				}
			}
		});
		waitPanel.startWaiting("<html>Please wait while " + taskName + "..." + "</html>");
		updateWindowTitleBlockingState(taskName);
		backgroundThread.start();
	}

	public void viewHelpFor(OperationDefinition definition) {
        
		String url = super.getHelpFor(definition);
	    
	    try {
			BrowserLauncher.openURL(url);
		} catch (Exception e) {
			reportException(e);
		}
	}

	public void viewHelp(String page) {
		try {
			super.getHelpUrl(page);
			BrowserLauncher.openURL(page);
		} catch (Exception e) {
			reportException(e);
		}
	}

	protected void garbageCollect() {
		System.gc();
		statusBar.updateMemoryIndicator();
	}

	public TaskManagerScreen getTaskManagerScreen() {
		TaskExecutor taskExecutor = Session.getSession().getServiceAccessor().getTaskExecutor();
		return new TaskManagerScreen(taskExecutor);
	}
	
	public void loadSession(boolean remote) {
		loadSession(remote, false, true);		
	}
	
	public void loadSession(boolean remote, boolean openExampleDir, boolean clear) {

		try {
			// create filechooser dialog
			final JFileChooser fileChooser = getSessionFileChooser(remote, openExampleDir);						

			int ret = fileChooser.showOpenDialog(this.getMainFrame());

			// user has selected a file
			if (ret == JFileChooser.APPROVE_OPTION) {
				File selectedFile = fileChooser.getSelectedFile();
				File sessionFile = null;
				String sessionId = null;
				String remoteSessionName = null;

				if (remote) {
					try {
						@SuppressWarnings("unchecked")
						List<DbSession> sessions = (List<DbSession>)fileChooser.getClientProperty("sessions");
						remoteSessionName = selectedFile.getPath().substring(ServerFile.SERVER_SESSION_ROOT_FOLDER.length()+1);
						DbSession session = getSessionManager().findSessionWithName(sessions, remoteSessionName);
						if (session == null) {
							showDialog("Session \"" + selectedFile + "\" not found", Severity.INFO, true);
							return;
						}
						sessionId = session.getDataId();
						if (sessionId == null) {
							// user didn't select anything
							showDialog("Session \"" + selectedFile + "\" not found", Severity.INFO, true);
							return;
						}

					} catch (Exception e) {
						reportException(e);
						throw new RuntimeException("internal error: URL or name from save dialog was invalid"); // should never happen
					}

				} else {
					// check that file exists
					if (!selectedFile.exists()) {
						DialogInfo info = new DialogInfo(Severity.INFO, "Could not open session file.", "File '" + selectedFile.getName() + "' not found.", "", Type.MESSAGE);
						ChipsterDialog.showDialog(this, info, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, true);
						return;
					}

					// check that the file is a session file
					if (!UserSession.isValidSessionFile(selectedFile)) {
						DialogInfo info = new DialogInfo(Severity.INFO, "Could not open session file.", "File '" + selectedFile.getName() + "' is not a valid session file.", "", Type.MESSAGE);
						ChipsterDialog.showDialog(this, info, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, true);
						return;
					}

					sessionFile = selectedFile;
				} 

				int xOffset = 0;

				try {
					// clear previous session 
					if (clear) {
						if (!clearSession()) {
							return; // loading cancelled
						}
					} else {
						xOffset = (int) getGraphPanel().getGraph().getGraphSize().getWidth();
					}

					// load the new session
					loadSession(sessionFile, sessionId, remote, false, false, xOffset);			
				} catch (MalformedURLException | FileBrokerException e) {
					reportException(e);
				}
			}
			menuBar.updateMenuStatus();
			
		} catch (Exception e) {
			if (remote) {
				DialogInfo info = new DialogInfo(Severity.ERROR, "Could not connect to server.", "Currently there is a problem in the network connection or the file server. During that time remote sessions are not accessible, but data can still be saved locally.", Exceptions.getStackTrace(e), Type.MESSAGE);
				ChipsterDialog.showDialog(this, info, DetailsVisibility.DETAILS_HIDDEN, true);
				return;
			} else {
				reportException(e);
			}
		}
	}
	
	public void loadSession(final File sessionFile, final String sessionId, final boolean isDataless, final boolean clearDeadTempDirs, final boolean isExampleSession, final Integer xOffset) {

		// start loading the session
		runBlockingTask("loading the session", new Runnable() {
			public void run() {
				getSessionManager().loadSessionAndWait(sessionFile, sessionId, isDataless, clearDeadTempDirs, isExampleSession, xOffset);
			}
		});
	}
	
	public void restoreSessionFrom(final File file) {
		
		runBlockingTask("loading the session", new Runnable() {
			public void run() {
				getSessionManager().restoreSessionAndWait(file);
			}
		});
	}
	
	public void saveSession(final SessionSavingMethod savingMethod) {

		// remote and local saves are quite different, first check which one this is
		final boolean remote;
		switch (savingMethod) {
		case INCLUDE_DATA_INTO_ZIP:
			remote = false;
			break;
		case UPLOAD_DATA_TO_SERVER:
			remote = true;
			break;
		default: 
			throw new IllegalArgumentException("internal error, not supported: " + savingMethod);
		}

		// create filechooser dialog
		try {
			JFileChooser fileChooser = getSessionFileChooser(remote, false);

			int ret = fileChooser.showSaveDialog(this.getMainFrame());

			// if was approved, then save it
			if (ret == JFileChooser.APPROVE_OPTION) {
				try {
					final File file;
					boolean exists;

					if (remote) {
						// use filename as it is (remote sessions use more human readable names)
						file = fileChooser.getSelectedFile();
						exists = false;

						@SuppressWarnings("unchecked")
						List<DbSession> sessions = (List<DbSession>)fileChooser.getClientProperty("sessions");
						for (DbSession session : sessions) {
							if (file.getName().equals(session.getName())) {
								exists = true;
								break;
							}						
						}

					} else {
						// add extension if needed
						file = fileChooser.getSelectedFile().getName().endsWith("." + UserSession.SESSION_FILE_EXTENSION) ? fileChooser.getSelectedFile() : new File(fileChooser.getSelectedFile().getCanonicalPath() + "." + UserSession.SESSION_FILE_EXTENSION);
						exists = file.exists();
					}

					// check if file (local or remote) exists
					if (exists) {
						int returnValue = JOptionPane.DEFAULT_OPTION;

						String message = "The file " + file.getName() + " already exists. Do you want to replace it?";

						Object[] options = { "Cancel", "Replace" };

						returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), message, "Confirm replace", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

						if (returnValue != 1) {
							return;
						}
					}

					// block GUI while saving
					runBlockingTask("saving session", new Runnable() {

						public void run() {

							// save
							getSessionManager().saveSessionAndWait(remote, file, file.getName());																				

							menuBar.updateMenuStatus();						
						}
					});
				} catch (Exception exp) {
					showErrorDialog("Saving session failed.", exp);
					return;
				}
			}
			menuBar.updateMenuStatus();

		} catch (Exception e) {
			if (remote) {
				DialogInfo info = new DialogInfo(Severity.ERROR, "Could not connect to server.", "Currently there is a problem in the network connection or the file server. During that time remote sessions are not accessible, but data can still be saved locally.", Exceptions.getStackTrace(e), Type.MESSAGE);
				ChipsterDialog.showDialog(this, info, DetailsVisibility.DETAILS_HIDDEN, true);
				return;
			} else {
				reportException(e);
			}
		}
	}	

	/**
	 * @return true if cleared, false if canceled
	 * @throws JMSException 
	 * @throws MalformedURLException 
	 */
	public boolean clearSession() throws MalformedURLException, FileBrokerException {

		if (!killUploadingTasks()) {
			return false;
		}
		
		int returnValue = JOptionPane.DEFAULT_OPTION;
		if (getSessionManager().hasUnsavedChanges()) {

			String message = "The current session contains unsaved changes.\nDo you want to clear it anyway?";
			Object[] options = { "Cancel", "Clear" };
			returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), message, "Clear session", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
		}

		if (!getSessionManager().hasUnsavedChanges() || returnValue == 1) {
			getSessionManager().clearSessionWithoutConfirming();			
			return true;
		}
		
		return false;
	}

	/**
	 * Method is called periodically to maintain state that cannot be maintained 
	 * in realtime. 
	 */
	public void checkFreeMemory() {
		statusBar.updateMemoryIndicator();
	}

	public VisualisationFrameManager getVisualisationFrameManager() {
		return visualisationFrameManager;
	}

	public void viewTasks() {
		statusBar.viewTasks();		
	}

	private void taskCountChanged(int newTaskCount, boolean attractAttention) {
		int completion = 0;
		if (newTaskCount > 0) {
			completion = taskExecutor.getTasks(true, false).iterator().next().getCompletionPercentage();
		}
		statusBar.taskCountChanged(newTaskCount, completion, attractAttention);		
	}

	public void refreshTaskList() {
		this.taskManagerScreen.refreshTasks();
	}

	public Screen getTaskListScreen() {
		return childScreens.get("TaskList");
	}

	public void showRenameView() {
		Variable renameVariable = new Variable(DataDetails.COMMAND, DataDetails.RENAME_COMMAND);
		
		setVisualisationMethod(
				VisualisationMethods.DATA_DETAILS, Arrays.asList(new Variable[] {renameVariable}), 
				getSelectionManager().getSelectedDataBeans(), 
				FrameType.MAIN); 
	}

	public void manageRemoteSessions() {
		
		final JFileChooser fileChooser = new RemoteSessionChooserFactory(this).getManagementChooser();							
		fileChooser.showOpenDialog(this.getMainFrame());
	}

	private void enableKeyboardShortcuts() {
		// add application wide keyboard shortcuts
		final HashMap<KeyStroke, Action> shortcutActionMap = new HashMap<KeyStroke, Action>();
					
		// find tool
		
		KeyStroke findKeyStroke = KeyStroke.getKeyStroke(KeyEvent.VK_F, Toolkit.getDefaultToolkit().getMenuShortcutKeyMask());
		shortcutActionMap.put(findKeyStroke,	
				new AbstractAction("FIND_TOOL") {
			@Override
			public void actionPerformed(ActionEvent e) {
				toolPanel.focusSearchField();
			}
		});
		// KeyStroke.toString() creates ugly "ctrl/meta pressed F", make it little bit more like
		// menu shortcuts
		String shortcut = findKeyStroke.toString().replace(" pressed ", "-");		
		toolPanel.setSearchFieldToolTipText("Find tool (" + shortcut + ")");
		
		// debug print
		
		shortcutActionMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_1, KeyEvent.CTRL_DOWN_MASK | KeyEvent.ALT_DOWN_MASK | KeyEvent.SHIFT_DOWN_MASK), 
				new AbstractAction("DEBUG_PRINT_SESSION") {
			@Override
			public void actionPerformed(ActionEvent e) {
				showDebugDialog(1);
			}
		});
		
		// open feedback session
		
		shortcutActionMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_O, 
				KeyEvent.CTRL_DOWN_MASK | KeyEvent.ALT_DOWN_MASK | 
				KeyEvent.SHIFT_DOWN_MASK),
				new AbstractAction("OPEN_LIGHTWEIGHT_SESSION") {
			@Override
			public void actionPerformed(ActionEvent e) {
				
			    String sessionId = JOptionPane.showInputDialog(
			    		mainFrame,
			    		"sessionId",
		                "Open feedback session",
		                JOptionPane.PLAIN_MESSAGE);		    
				
			    if (sessionId != null) {
			    	loadSession(null, sessionId, true, false, false, null);
			    }
			}
		});
		KeyboardFocusManager kfm = 
				KeyboardFocusManager.getCurrentKeyboardFocusManager();
		kfm.addKeyEventDispatcher(new KeyEventDispatcher() {
			@Override
			public boolean dispatchKeyEvent(KeyEvent e) {
				KeyStroke keyStroke = KeyStroke.getKeyStrokeForEvent(e);
				if (shortcutActionMap.containsKey(keyStroke) ) {
					final Action a = shortcutActionMap.get(keyStroke);
					final ActionEvent ae = new ActionEvent(e.getSource(), e.getID(), 
							null );
					SwingUtilities.invokeLater( new Runnable() {
						@Override
						public void run() {
							a.actionPerformed(ae);
						}
					});
					return true;
				}
				return false;
			}
		});
	}
	
	public GraphPanel getGraphPanel() {
		return graphPanel;
	}
}
