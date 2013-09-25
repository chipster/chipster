package fi.csc.microarray.client;

import java.awt.BorderLayout;
import java.awt.CardLayout;
import java.awt.Component;
import java.awt.KeyEventDispatcher;
import java.awt.KeyboardFocusManager;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;
import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.BorderFactory;
import javax.swing.Icon;
import javax.swing.JComponent;
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
import fi.csc.microarray.client.dataimport.ImportUtils.FileLoaderProcess;
import fi.csc.microarray.client.dataimport.table.InformationDialog;
import fi.csc.microarray.client.dataview.DetailsPanel;
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
import fi.csc.microarray.client.dialog.SnapshotAccessory;
import fi.csc.microarray.client.dialog.URLImportDialog;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.ToolPanel;
import fi.csc.microarray.client.screen.ChildScreenPool;
import fi.csc.microarray.client.screen.HistoryScreen;
import fi.csc.microarray.client.screen.Screen;
import fi.csc.microarray.client.screen.ShowSourceScreen;
import fi.csc.microarray.client.screen.TaskManagerScreen;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.session.SessionLoader.LoadMethod;
import fi.csc.microarray.client.session.UserSession;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.waiting.WaitGlassPane;
import fi.csc.microarray.client.workflow.WorkflowManager;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.ContentType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.Traversal;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.exception.ErrorReportAsException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;
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

	private static final int METADATA_FETCH_TIMEOUT_SECONDS = 15;
	private static final long SLOW_VISUALISATION_LIMIT = 5 * 1000;
	private static final long VERY_SLOW_VISUALISATION_LIMIT = 20 * 1000;

	/**
	 * Logger for this class
	 */
	private static Logger logger;

	private JFrame mainFrame = null;
	private JPanel rightSideViewChanger = null;
	private JPanel leftSideContentPane = null;
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
	private SimpleInternalFrame detailsFrame;

	private ChildScreenPool childScreens;
	private TreePanel tree;
	private DetailsPanel details;
	private GraphPanel graphPanel;
	private ToolPanel toolPanel;
	private VisualisationFrameManager visualisationFrameManager;
	private HistoryScreen historyScreen;

	private SplashScreen splashScreen;
	private ClientListener clientListener;
	private WaitGlassPane waitPanel = new WaitGlassPane();
	
	private static float fontSize = VisualConstants.DEFAULT_FONT_SIZE;

	private JFileChooser importExportFileChooser;
	private JFileChooser sessionFileChooser;
	private JFileChooser workflowFileChooser;

	public SwingClientApplication(ClientListener clientListener, AuthenticationRequestListener overridingARL, String module, boolean isStandalone)
	        throws MicroarrayException, IOException, IllegalConfigurationException {

		super(isStandalone, overridingARL);
		
		// this had to be delayed as logging is not available before loading configuration
		logger = Logger.getLogger(SwingClientApplication.class);
		
		if (!SwingUtilities.isEventDispatchThread()) {
			logger.error(new MicroarrayException("SwingClientApplication was created outside the Event Dispatch Thread."));
			System.exit(1);
		}
		
		this.clientListener = clientListener;

        // set the module that user wants to load
        this.requestedModule = module;

        // show splash screen
		splashScreen = new SplashScreen(VisualConstants.SPLASH_SCREEN);
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
					initialiseApplication();

				} catch (Exception e) {
					reportInitalisationErrorThreadSafely(e);
				}
			}
		});
		t.start();
	}
	
	private void reportInitalisationErrorThreadSafely(final Exception e) {
		SwingUtilities.invokeLater(new Runnable() {
			
			@Override
			public void run() {

				showDialog("Starting Chipster failed.", "There could be a problem with the network connection, or the remote services could be down. " +
						"Please see the details below for more information about the problem.\n\n" + 
						"Chipster also fails to start if there has been a version update with a change in configurations. In such case please delete Chipster application settings directory.",
						Exceptions.getStackTrace(e), Severity.ERROR, false, ChipsterDialog.DetailsVisibility.DETAILS_HIDDEN,
						new PluginButton() {
					@Override
					public void actionPerformed() {
						try {
							new SwingClientApplication(getShutdownListener(), null, null, true);
							
						} catch (Exception e) {
							// ignore
						}
					}
					@Override
					public String getText() {
						return "Start standalone";
					}
				});
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

	protected void initialiseGUIThreadSafely(final File mostRecentDeadTempDirectory) {
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
			throw new MicroarrayException("metadata was not received (analyser not functional?)");
		}

		// initialize the main frame
		this.mainFrame = new JFrame();
		updateWindowTitle();
		childScreens = new ChildScreenPool(mainFrame);
		Frames frames = new Frames(mainFrame);
		Session.getSession().setFrames(frames);

		// Sets look 'n' feel
		setPlastic3DLookAndFeel(mainFrame);

		// set location
		mainFrame.setLocationByPlatform(true);

		// initialise joblist popup menu
		// do this method before getStatusBar to avoid null pointer exception
		this.taskManagerScreen = this.getTaskManagerScreen();

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

		// initialise left side
		leftSideContentPane = new JPanel(new BorderLayout());
		details = new DetailsPanel(leftSideContentPane);
		leftSideContentPane.setBorder(BorderFactory.createEmptyBorder());

		/* Initialize tree and graph */
		
		//moved to getTreeFrame
//		this.tree = new TreePanel(manager.getRootFolder());
		this.graphPanel = new GraphPanel();

		treeFrame = getTreeFrame();
		graphFrame = getGraphFrame();

		leftSplit = new JSplitPane(JSplitPane.VERTICAL_SPLIT, treeFrame, graphFrame);
		leftSplit.setDividerLocation(VisualConstants.TREE_PANEL_HEIGHT);
		leftSplit.setResizeWeight(0.1);

		detailsFrame = getDetailsFrame();

		leftSideContentPane.add(leftSplit, BorderLayout.CENTER);
		leftSideContentPane.add(detailsFrame, BorderLayout.SOUTH);

		rightSideViewChanger = new JPanel(new BorderLayout());

		if (this.isStandalone()) {
			rightSideViewChanger.add(visualisationArea, BorderLayout.CENTER);
		} else {	
			rightSideViewChanger.add(rightSplit, BorderLayout.CENTER);
		}
		rightSideViewChanger.setBorder(BorderFactory.createEmptyBorder());

		// construct the whole main content pane
		mainSplit = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, leftSideContentPane, rightSideViewChanger);
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
		mainFrame.setIconImage(VisualConstants.APPLICATION_ICON.getImage());
		mainFrame.pack();
		mainFrame.setExtendedState(JFrame.MAXIMIZED_BOTH);
		mainFrame.setVisible(true);
		mainFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		
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
		customiseFocusTraversal();
		restoreDefaultView();
		enableKeyboardShortcuts();
		
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

	private void customiseFocusTraversal() throws MicroarrayException {
		Vector<Component> order = new Vector<Component>();
		order.addAll(tree.getFocusComponents());
		order.addAll(toolPanel.getFocusComponents());
		order.addAll(visualisationFrameManager.getFocusComponents());

		getMainFrame().setFocusTraversalPolicy(new ClientFocusTraversalPolicy(order));

	}

	
	private String windowTitleJobPrefix = null;
	private String windowTitleBlockingPrefix = null;

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

	public SimpleInternalFrame getDetailsFrame() throws MicroarrayException {
		if (detailsFrame == null) {

			class DetailsFrame extends SimpleInternalFrame implements PropertyChangeListener {
				public DetailsFrame() {
					super("Notes for dataset");
				}

				public void propertyChange(PropertyChangeEvent evt) {
					if (evt instanceof DatasetChoiceEvent) {
						DatasetChoiceEvent dce = (DatasetChoiceEvent) evt;
						if (dce.getNewValue() != null) {
							this.setTitle("Notes for dataset " + dce.getNewValue());
						} else {
							this.setTitle("Notes for dataset");
						}
					}
				}
			}

			DetailsFrame detailsFrameWithListener = new DetailsFrame();
			addClientEventListener(detailsFrameWithListener);

			this.detailsFrame = detailsFrameWithListener;

			detailsFrame.add(details);
			return detailsFrame;

		} else {
			return detailsFrame;
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

	@Override
	public void importGroup(final Collection<ImportItem> datas, final String folderName) {

		runBlockingTask("importing files", new Runnable() {

			public void run() {
				DataBean lastGroupMember = null;

				try {

					for (ImportItem item : datas) {

						String dataSetName = item.getOutput().getName();
						ContentType contentType = item.getType();
						Object dataSource = item.getInput();


						// Selects folder where data is imported to, or creates a
						// new one
						DataFolder folder = initializeFolderForImport(folderName);

						// get the InputStream for the data source
						InputStream input;

						
						// create the DataBean
						DataBean data = null;

						if (dataSource instanceof File) {
							data = manager.createDataBean(dataSetName, (File) dataSource);
							
						} else if (dataSource instanceof URL) {
							// TODO Not used anymore, URL-files are saved to the
							// temp file
							URL url = (URL) dataSource;
							try {
								input = url.openStream();
								manager.createDataBean(dataSetName, input);
							} catch (FileNotFoundException fnfe) {
								SwingUtilities.invokeAndWait(new Runnable() {
									public void run() {
										showDialog("File not found.", null, "File not found. Check that the typed URL is pointing to a valid location", Severity.ERROR, false);
									}
								});
								break;

							} catch (IOException ioe) {
								SwingUtilities.invokeAndWait(new Runnable() {
									public void run() {
										showDialog("Import failed.", null, "Error occured while importing data from URL", Severity.ERROR, false);
									}
								});
								break;
							}

						} else if (dataSource instanceof InputStream) {
							logger.info("loading data from a plain stream, caching can not be used!");
							input = (InputStream) dataSource;

						} else {
							throw new IllegalArgumentException("unknown dataSource type: " + dataSource.getClass().getSimpleName());
						}

						// make sure that the new bean is not null
						if (data == null) {
							SwingUtilities.invokeAndWait(new Runnable() {
								public void run() {
									showDialog("Importing dataset filed.", null, "Created DataBean was null.", Severity.WARNING, false);
								}
							});
							return;
						}
						
						// set the content type
						data.setContentType(contentType);

						// add the operation (all databeans have their own import
						// operation
						// instance, it would be nice if they would be grouped)
						Operation importOperation = new Operation(OperationDefinition.IMPORT_DEFINITION, new DataBean[] { data });
						data.setOperationRecord(new OperationRecord(importOperation));

						// data is ready now, make it visible
						folder.addChild(data);

						// Create group links only if both datas are raw type
						if (lastGroupMember != null && ChipsterInputTypes.hasRawType(lastGroupMember) && ChipsterInputTypes.hasRawType(data)) {

							DataBean targetData = data;

							// Link new data to all group linked datas of given cell
							for (DataBean sourceData : lastGroupMember.traverseLinks(new Link[] { Link.GROUPING }, Traversal.BIDIRECTIONAL)) {
								logger.debug("Created GROUPING link between " + sourceData.getName() + " and " + targetData.getName());
								createLink(sourceData, targetData, DataBean.Link.GROUPING);
							}

							// Create link to the given cell after looping to avoid
							// link duplication
							createLink(lastGroupMember, targetData, DataBean.Link.GROUPING);
						}

						lastGroupMember = data;

					}

				} catch (Exception e) {
					throw new RuntimeException(e);
				}

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

	
	
	
	@Override
	public File saveWorkflow() {

		try {
			JFileChooser fileChooser = this.getWorkflowFileChooser();
			int ret = fileChooser.showSaveDialog(this.getMainFrame());
			if (ret == JFileChooser.APPROVE_OPTION) {
				File selected = fileChooser.getSelectedFile();
				File newFile = selected.getName().endsWith(WorkflowManager.SCRIPT_EXTENSION) ? selected : new File(selected.getCanonicalPath() + "." + WorkflowManager.SCRIPT_EXTENSION);

				workflowManager.saveSelectedWorkflow(newFile);
				unsavedChanges = false;
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

	public void runWorkflow(URL workflowScript) {
		runWorkflow(workflowScript, null);
	}

	public void runWorkflow(URL workflowScript, final AtEndListener atEndListener) {
		workflowManager.runScript(workflowScript, atEndListener);
	}


	@Override
	public File openWorkflow() {

		try {
			JFileChooser fileChooser = this.getWorkflowFileChooser();
			int ret = fileChooser.showOpenDialog(this.getMainFrame());
			if (ret == JFileChooser.APPROVE_OPTION) {
				runWorkflow(fileChooser.getSelectedFile().toURI().toURL());

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

	public void showDetailsFor(DataBean data) {
		details.setViewedData(data);
	}

	public Icon getIconFor(DataItem element) {
		if (element instanceof DataFolder) {
			return VisualConstants.ICON_TYPE_FOLDER;
		} else {
			DataBean bean = (DataBean) element;
			if (bean.queryFeatures("/phenodata").exists()) {
				return VisualConstants.ICON_TYPE_PHENODATA;
			} else {
				return bean.getContentType().getIcon();
			}
		}
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
		
		//
		boolean userFixable = task.getState() == State.FAILED_USER_ERROR && task.getErrorMessage() != null && !task.getErrorMessage().equals("");
		
		// user-friendly message
		if (userFixable) {
			title = task.getErrorMessage();
			message = task.getNamePrettyPrinted() + " was stopped. ";
		} 
		
		// generic message
		else {
			title = task.getNamePrettyPrinted() + " did not complete successfully. ";
			message = "You may have used a tool or parameters which are unsuitable for the selected dataset, or " + "there might be a bug in the analysis tool itself.\n\n" + "The details below may provide hints about the problem.";
		}		

		// details
		String details = "";
		if (task.getErrorMessage() != null) {
			details += task.getErrorMessage() + "\n\n";
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
		        "If problem persist, please check that your data is valid. For more information open the details panel below.", null);
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

			if (ErrorDialogUtils.getMessage(me.getExtraInfo()) != null) {
				details += ErrorDialogUtils.getMessage(me.getExtraInfo()) + "\n";
			}

			if (ErrorDialogUtils.getScreenOutput(me.getExtraInfo()) != null) {
				details += ErrorDialogUtils.getScreenOutput(me.getExtraInfo());
			}

			dialogInfo.setDetails(details);

		} else {
			// use stack trace as details
			dialogInfo.setDetails(Exceptions.getStackTrace(e));
		}

		// show dialog
		ChipsterDialog.showDialog(this, dialogInfo, ChipsterDialog.DetailsVisibility.DETAILS_HIDDEN, false);

		// we'll always output these to console and log for traceability and
		// easier IDE navigation
		e.printStackTrace();
		if (logger != null) {
			logger.error(Exceptions.getStackTrace(e));
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
		return getDefaultVisualisationForSelection() != VisualisationMethod.NONE;
	}

	/**
	 * Gets default visualisation method for selected databeans. The method is
	 * selected by following steps:
	 * 
	 * <ol>
	 * <li>If no dataset is selected, return
	 * <code>VisualisationMethod.NONE</code> </li>
	 * <li>If only one dataset is selected, return the default method for the
	 * data </li>
	 * </li>
	 * <li>If multiple datasets are selected, check the best method for each
	 * dataset. If the best method is same for all selected datasets and it can
	 * be used with multiple data, the best method is returned. </li>
	 * <li>If the best method is not same for all of the datas, try to find
	 * just some method which is suitable for all datas and can be used with
	 * multiple datasets. </li>
	 * <li>If there were no method to fill the requirements above, return
	 * <code>VisualisationMethod.NONE</code> </li>
	 * 
	 * @return default visualisation method which is suitable for all selected
	 *         datasets
	 */
	private VisualisationMethod getDefaultVisualisationForSelection() {
		logger.debug("getting default visualisation");
		if (getSelectionManager().getSelectedDataBeans() == null || getSelectionManager().getSelectedDataBeans().size() == 0) {
			return VisualisationMethod.NONE;
		}

		try {
			List<DataBean> beans = getSelectionManager().getSelectedDataBeans();

			if (beans.size() == 1) {
				return Session.getSession().getVisualisations().getDefaultVisualisationFor(beans.get(0));
			} else if (beans.size() > 1)
				for (VisualisationMethod method : Session.getSession().getVisualisations().getOrderedDefaultCandidates()) {
					if (method == VisualisationMethod.NONE || !method.getHeadlessVisualiser().isForMultipleDatas()) {
						continue;
					}
					if (method.isApplicableTo(beans)) {
						return method;
					}
				}

			/*
			 * 
			 * VisualisationMethod defaultMethodForDatas = null; // First, try
			 * to find best suitable visualisation for all for (DataBean bean :
			 * beans) { VisualisationMethod method = new
			 * BioBean(bean).getDefaultVisualisation(); if
			 * (defaultMethodForDatas == null &&
			 * VisualisationMethod.isApplicableForMultipleDatas(method)) {
			 * defaultMethodForDatas = method; } else { if
			 * (defaultMethodForDatas != method) { // Searching for best method
			 * for all failed defaultMethodForDatas = null; logger.debug("Method " +
			 * method + " can not be used to visualise selected datas"); break; } } }
			 * 
			 * if (defaultMethodForDatas != null) { // Visualise datas if the
			 * best method was found logger.debug("Method " +
			 * defaultMethodForDatas + " will be used to visualise selected
			 * datas"); return defaultMethodForDatas; } // Keep looking for
			 * suitable visualisation DataBean firstData = beans.get(0);
			 * 
			 * for (VisualisationMethod method :
			 * VisualisationMethod.getApplicableForMultipleDatas()) { if (method ==
			 * VisualisationMethod.NONE) { continue; }
			 * 
			 * if (method.isApplicableTo(firstData)) { // The method is
			 * applicable to one of the selected datas // Check that the same
			 * method is applicable to the other // datasets too boolean
			 * isSuitableMethod = true; for (DataBean otherData : beans) { if
			 * (otherData.equals(firstData)) { continue; }
			 * 
			 * if (!method.isApplicableTo(otherData)) { isSuitableMethod =
			 * false; logger.debug("Method " + method + " can not be used to
			 * visualise selected datas"); break; } }
			 * 
			 * if (isSuitableMethod) { logger.debug("Method " + method + " will
			 * be used to visualise selected datas"); return method; } } }
			 */
			return VisualisationMethod.NONE;

		} catch (Exception e) {
			reportException(e);
			return VisualisationMethod.NONE;
		}
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
		deleteDatasWithoutConfirming(datas);
	}
	
	public void deleteDatasWithoutConfirming(DataItem... datas) {
	
		// check that we have something to delete
		if (datas.length == 0) {
			return; // no selection, do nothing
		}		
		
		// remove all selections
		getSelectionManager().clearAll(true, this);

		// do actual delete
		for (DataItem data : datas) {
			manager.delete(data);
		}

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
			List<File> files = new ArrayList<File>();

			for (File file : fc.getSelectedFiles()) {
				files.add(file);
			}

			ImportSession importSession = new ImportSession(ImportSession.Source.FILES, files, access.getImportFolder(), access.skipActionChooser());
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

			File file = ImportUtils.createTempFile(ImportUtils.URLToFilename(selectedURL), ImportUtils.getExtension(ImportUtils.URLToFilename(selectedURL)));

			ImportUtils.getURLFileLoader().loadFileFromURL(selectedURL, file, importFolder, urlImportDlg.isSkipSelected());
		}
	}
	
	protected void quit() {
		int returnValue = JOptionPane.DEFAULT_OPTION;

		// Check the running tasks
		if (taskExecutor.getRunningTaskCount() > 0) {
			String message = "";
			if (taskExecutor.getRunningTaskCount() == 1) {
				message += "There is a running task.  Are you sure you want to cancel the running task?";
			} else {
				message += "There are " + taskExecutor.getRunningTaskCount() + " running tasks. " + "Are you sure you want to cancel all running tasks?";
			}

			Object[] options = { "Cancel running tasks", "Cancel" };

			returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), message, "Confirm close", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

			if (returnValue == JOptionPane.YES_OPTION) {
				taskExecutor.killAll();
			} else {
				return;
			}
		}

		// Check for unsaved changes
		returnValue = JOptionPane.DEFAULT_OPTION;

		if (unsavedChanges) {

			Object[] options = { "Save and close", "Close without saving", "Cancel" };

			returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), "Do you want the session to be saved before closing Chipster?", "Confirm close", JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

			if (returnValue == 0) {
				try {
					saveSessionAndQuit();
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
	
	public void quitImmediately() {

		// hide immediately to look more reactive...
		childScreens.disposeAll();
		mainFrame.setVisible(false);

		super.quit();

		// this closes the application
		mainFrame.dispose();
		System.exit(0);
	}

	
	public static void startStandalone(String module) throws IOException {
		try {
			DirectoryLayout.initialiseStandaloneClientLayout();
			Configuration config = DirectoryLayout.getInstance().getConfiguration();
			config.getRootModule().getModule("messaging").getEntry("broker-host").setValue("(none)");
			config.getRootModule().getModule("messaging").getEntry("broker-protocol").setValue("");
			config.getRootModule().getModule("messaging").getEntry("broker-port").setValue("0");
			config.getRootModule().getModule("security").getEntry("username").setValue("");
			config.getRootModule().getModule("security").getEntry("password").setValue("");
					
		} catch (IllegalConfigurationException e) {
			reportIllegalConfigurationException(e);
		}

		ClientListener shutdownListener = getShutdownListener();
		
		try {
			new SwingClientApplication(shutdownListener, null, module, true);
			
		} catch (Throwable t) {
			t.printStackTrace();
			if (logger != null) {
				logger.error(Exceptions.getStackTrace(t));
			}
		}

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
			new SwingClientApplication(shutdownListener, null, module, false);
			
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

	@Override
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
	public void setMaximisedVisualisationMode(boolean maximisedVisualisationMode) {
		showMaximisedVisualisation(maximisedVisualisationMode);
	}

	@Override
	public void setVisualisationMethod(VisualisationMethod method, List<Variable> variables, List<DataBean> datas, FrameType target) {

		if (method == null || datas == null) {
			super.setVisualisationMethod(VisualisationMethod.NONE, null, null, target);
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

	private JFileChooser getSessionFileChooser(JComponent accessory) {
		if (sessionFileChooser == null) {
			sessionFileChooser = ImportUtils.getFixedFileChooser();

			String[] extensions = { UserSession.SESSION_FILE_EXTENSION };
			sessionFileChooser.setFileFilter(new GeneralFileFilter("Chipster Session (." + UserSession.SESSION_FILE_EXTENSION + ")", extensions));
			sessionFileChooser.setSelectedFile(new File("session." + UserSession.SESSION_FILE_EXTENSION));
			sessionFileChooser.setAcceptAllFileFilterUsed(false);
			sessionFileChooser.setMultiSelectionEnabled(false);

		}
		sessionFileChooser.setAccessory(accessory);

		fixFileChooserFontSize(sessionFileChooser);

		return sessionFileChooser;
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
	 * 
	 * @param areaRequestHandler
	 */
	public void openImportTool(ImportSession importSession) {
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
        String url = definition.getHelpURL();
	    if (url != null && !url.isEmpty()) {
	        // Link is stored in operation definition
	        url = definition.getHelpURL();
	    } else {
	        // Mostly for microarray
	        // TODO: consider refactoring so that url is stored in definition
	        // and this "else" branch is not needed
	        url = HelpMapping.mapToHelppage(definition);
	    }
	    
	    try {
			BrowserLauncher.openURL(url);
		} catch (Exception e) {
			reportException(e);
		}
	}

	public void viewHelp(String page) {
		try {
			if (!page.startsWith(HelpMapping.MANUAL_ROOT)) {
				page = HelpMapping.MANUAL_ROOT + page;
			}
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

	public void createLink(DataBean source, DataBean target, Link type) {
		source.addLink(type, target);
	}

	public void removeLink(DataBean source, DataBean target, Link type) {
		source.removeLink(type, target);
	}

	@Override
	public void showImportToolFor(File file, String destinationFolder, boolean skipActionChooser) {
		ImportSession importSession = new ImportSession(ImportSession.Source.FILES, new File[] { file }, destinationFolder, skipActionChooser);
		openImportTool(importSession);
	}


	@Override
	public void loadSessionFrom(File file) {
		loadSessionImpl(file, LoadMethod.NORMAL);
	}

	@Override
	public void restoreSessionFrom(File file) {
		loadSessionImpl(file, LoadMethod.LIGHTWEIGHT);
	}

	@Override
	public void loadSessionFrom(URL url) {
		try {
			final File tempFile = ImportUtils.createTempFile(ImportUtils.URLToFilename(url), ImportUtils.getExtension(ImportUtils.URLToFilename(url)));
			InformationDialog info = new InformationDialog("Loading session", "Loading session from the specified URL", null);
			
			FileLoaderProcess fileLoaderProcess = new FileLoaderProcess(tempFile, url, info) {
				@Override
				protected void postProcess() {
					loadSessionImpl(tempFile, LoadMethod.NORMAL);
				};
			};			
			fileLoaderProcess.runProcess();
			
		} catch (IOException e) {
			reportException(e);
		}
	}
	
	
	@Override
	public void loadSession() {
		loadSession(LoadMethod.NORMAL);
	}
	
	public void loadSession(LoadMethod loadMethod) {

		SnapshotAccessory accessory = new SnapshotAccessory();
		final JFileChooser fileChooser = getSessionFileChooser(accessory);
		int ret = fileChooser.showOpenDialog(this.getMainFrame());
		
		// user has selected a file
		if (ret == JFileChooser.APPROVE_OPTION) {
			File sessionFile = fileChooser.getSelectedFile();
			
			// check that file exists
			if (!sessionFile.exists()) {
				DialogInfo info = new DialogInfo(Severity.INFO, "Could not open session file.", "File '" + sessionFile.getName() + "' not found.", "", Type.MESSAGE);
				ChipsterDialog.showDialog(this, info, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, true);
				return;
			}
			
			// check that the file is a session file
			if (!UserSession.isValidSessionFile(sessionFile)) {
				DialogInfo info = new DialogInfo(Severity.INFO, "Could not open session file.", "File '" + sessionFile.getName() + "' is not a valid session file.", "", Type.MESSAGE);
				ChipsterDialog.showDialog(this, info, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, true);
				return;
			}
			
			// clear previous session 
			if (accessory.clearSession()) {
				if (!clearSession()) {
					return; // loading cancelled
				}
			}								

			// load the new session
			loadSessionImpl(fileChooser.getSelectedFile(), loadMethod);		
		}
		menuBar.updateMenuStatus();
	}

	private void loadSessionImpl(final File sessionFile, final LoadMethod loadMethod) {
		
		// check that it's a valid session file 
		if (!UserSession.isValidSessionFile(sessionFile)) {
			DialogInfo dialogInfo = new DialogInfo(Severity.INFO, "Could not open session file.", "The given file is not a valid session file.", "");
			ChipsterDialog.showDialog(this, dialogInfo, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, true);
			return;
		}
		
		// start loading the session
		runBlockingTask("loading the session", new Runnable() {
			public void run() {						
					
				/* If there wasn't data or it was just cleared, there is no need to warn about
				 * saving after opening session. However, if there was datasets already, combination
				 * of them and new session can be necessary to save. This has to set after the import, because 
				 */
				boolean somethingToSave = manager.databeans().size() != 0;

				manager.loadSession(sessionFile, loadMethod);

				unsavedChanges = somethingToSave;
				
				// should also be done for feedback sessions?
				if (loadMethod.equals(LoadMethod.LIGHTWEIGHT)) {
					clearDeadTempDirectories();
				}
			}
		});
	}
	
	
	@Override
	public void saveSession() {
		saveSession(false);
	}

	public void saveSessionAndQuit() {
		saveSession(true);
	}
	
	public void saveSession(final boolean quit) {

		JFileChooser fileChooser = getSessionFileChooser(null);
		int ret = fileChooser.showSaveDialog(this.getMainFrame());
		
		if (ret == JFileChooser.APPROVE_OPTION) {
			try {
				final File file = fileChooser.getSelectedFile().getName().endsWith("." + UserSession.SESSION_FILE_EXTENSION) ? fileChooser.getSelectedFile() : new File(fileChooser.getSelectedFile().getCanonicalPath() + "." + UserSession.SESSION_FILE_EXTENSION);

				if (file.exists()) {
					int returnValue = JOptionPane.DEFAULT_OPTION;

					String message = "The file " + file.getCanonicalPath() + " already exists. Do you want " + "to replace it?";

					Object[] options = { "Cancel", "Replace" };

					returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), message, "Confirm replace", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

					if (returnValue != 1) {
						return;
					}

				}

				// save
				runBlockingTask("saving session", new Runnable() {

					public void run() {

						// save
						boolean saveSuccessful;
						saveSuccessful = getDataManager().saveSession(file);

						if (saveSuccessful) {

							// quit
							if (quit) {
								quitImmediately();
							}

							menuBar.updateMenuStatus();
							unsavedChanges = false;
						}						
					}
				});
			} catch (Exception exp) {
				showErrorDialog("Saving session failed.", exp);
				return;
			}
		}
		menuBar.updateMenuStatus();
	}

	/**
	 * @return true if cleared, false if canceled
	 */
	public boolean clearSession() {

		int returnValue = JOptionPane.DEFAULT_OPTION;
		if (unsavedChanges) {

			String message = "The current session contains unsaved changes.\nDo you want to clear it anyway?";
			Object[] options = { "Cancel", "Clear" };
			returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), message, "Clear session", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
		}

		if (!unsavedChanges || returnValue == 1) {
			this.deleteDatasWithoutConfirming(manager.getRootFolder());
			unsavedChanges = false;
			return true;
		}
		return false;
	}

	@Override
	public void checkFreeMemory() {
		statusBar.updateMemoryIndicator();
	}

	public DataManager getDataManager() {
		return manager;
	}

	public VisualisationFrameManager getVisualisationFrameManager() {
		return visualisationFrameManager;
	}

	@Override
	public void flipTaskListVisibility(boolean closeIfVisible) {
		statusBar.flipTaskListVisibility(closeIfVisible);
		
	}

	@Override
	protected void taskCountChanged(int newTaskCount, boolean attractAttention) {
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
	
	private void enableKeyboardShortcuts() {
		// add application wide keyboard shortcuts
		final HashMap<KeyStroke, Action> shortcutActionMap = new 
				HashMap<KeyStroke, Action>();
		shortcutActionMap.put(KeyStroke.getKeyStroke(KeyEvent.VK_O, 
				KeyEvent.CTRL_DOWN_MASK | KeyEvent.ALT_DOWN_MASK | 
				KeyEvent.SHIFT_DOWN_MASK),
				new AbstractAction("OPEN_LIGHTWEIGHT_SESSION") {
			@Override
			public void actionPerformed(ActionEvent e) {
				loadSession(LoadMethod.FEEDBACK);
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
	
	
}
