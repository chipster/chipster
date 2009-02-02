package fi.csc.microarray.client;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;
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

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.ErrorReportAsException;
import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.dataimport.ImportScreen;
import fi.csc.microarray.client.dataimport.ImportSession;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.dataview.DetailsPanel;
import fi.csc.microarray.client.dataview.GraphPanel;
import fi.csc.microarray.client.dataview.TreePanel;
import fi.csc.microarray.client.dialog.ChipsterDialog;
import fi.csc.microarray.client.dialog.ClipboardImportDialog;
import fi.csc.microarray.client.dialog.DialogInfo;
import fi.csc.microarray.client.dialog.ErrorDialogUtils;
import fi.csc.microarray.client.dialog.ImportSettingsAccessory;
import fi.csc.microarray.client.dialog.SnapshotAccessory;
import fi.csc.microarray.client.dialog.URLImportDialog;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationPanel;
import fi.csc.microarray.client.screen.ChildScreenPool;
import fi.csc.microarray.client.screen.HistoryScreen;
import fi.csc.microarray.client.screen.Screen;
import fi.csc.microarray.client.screen.ShowSourceScreen;
import fi.csc.microarray.client.screen.TaskManagerScreen;
import fi.csc.microarray.client.selection.DatasetChoiceEvent;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.waiting.WaitGlassPane;
import fi.csc.microarray.client.workflow.WorkflowManager;
import fi.csc.microarray.databeans.ContentType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.Traversal;
import fi.csc.microarray.databeans.fs.FSSnapshottingSession;
import fi.csc.microarray.description.VVSADLParser.ParseException;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.auth.ClientLoginListener;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;
import fi.csc.microarray.util.BrowserLauncher;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.GeneralFileFilter;
import fi.csc.microarray.util.SplashScreen;
import fi.csc.microarray.util.config.ConfigurationLoader.OldConfigurationFormatException;

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

	/**W
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
	private OperationPanel operationsPanel;
	private VisualisationFrameManager visualisationFrameManager;
	private HistoryScreen historyScreen;

	private SplashScreen splashScreen;
	private ClientListener clientListener;
	private AuthenticationRequestListener overridingARL;
	private WaitGlassPane waitPanel = new WaitGlassPane();
	
	private static float fontSize = VisualConstants.DEFAULT_FONT_SIZE;

	private boolean unsavedChanges = false;

	private JFileChooser importExportFileChooser;
	private JFileChooser snapshotFileChooser;
	private JFileChooser workflowFileChooser;

	public SwingClientApplication(ClientListener clientListener, AuthenticationRequestListener overridingARL) throws MicroarrayException, IOException {

		this.clientListener = clientListener;
		this.overridingARL = overridingARL;

		splashScreen = new SplashScreen(VisualConstants.SPLASH_SCREEN);
		reportInitialisation("Initialising " + ApplicationConstants.APPLICATION_TITLE, true);

		// we want to close the splash screen exception occurs
		try {
			initialiseApplication();
		} catch (Exception e) {
			splashScreen.close();
			throw new MicroarrayException(e);
		}

		// this had to be delayed as logging is not available before loading
		// configuration
		logger = Logger.getLogger(SwingClientApplication.class);
	}

	public void reportInitialisation(String report, boolean newline) {
		if (newline) {
			splashScreen.writeLine(report);
		} else {
			splashScreen.write(report);
		}

	}

	protected void initialiseGUI() throws MicroarrayException, IOException {

		// assert state of initialisation
		try {
			definitionsInitialisedLatch.await(METADATA_FETCH_TIMEOUT_SECONDS, TimeUnit.SECONDS);
		} catch (InterruptedException e) {
			throw new RuntimeException(e);
		}
		if (parsedCategories == null) {
			throw new MicroarrayException("metadata was not received (analyser not functional?)");
		}

		// initialize the main frame
		mainFrame = new JFrame();
		updateWindowTitle(0); // zero jobs at start
		childScreens = new ChildScreenPool(mainFrame);

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
			operationsPanel = new OperationPanel(parsedCategories);
		} catch (ParseException e) {
			logger.error("VVSADL parse failed", e);
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
		this.tree = new TreePanel(manager.getRootFolder());
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
		rightSideViewChanger.add(rightSplit, BorderLayout.CENTER);
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

		// Remember changes to confirm close only when necessary
		manager.addDataChangeListener(new DataChangeListener() {
			public void dataChanged(DataChangeEvent event) {
				unsavedChanges = true;
			}
		});

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

		customiseFocusTraversal();
	}

	private void customiseFocusTraversal() throws MicroarrayException {
		Vector<Component> order = new Vector<Component>();
		order.addAll(tree.getFocusComponents());
		order.addAll(operationsPanel.getFocusComponents());
		order.addAll(visualisationFrameManager.getFocusComponents());

		getMainFrame().setFocusTraversalPolicy(new ClientFocusTraversalPolicy(order));

	}

	public void updateWindowTitle(Integer jobCount) {
		String jobString = jobCount > 0 ? jobCount + " jobs / " : "";
		this.mainFrame.setTitle(jobString + ApplicationConstants.APPLICATION_TITLE);
	}

	public SimpleInternalFrame getOperationsFrame() {
		if (operationsFrame == null) {
			operationsFrame = new SimpleInternalFrame("Analysis tools");
			operationsFrame.setContent(operationsPanel);
		}
		return operationsFrame;
	}

	public SimpleInternalFrame getGraphFrame() throws MicroarrayException {
		if (graphFrame == null) {
			graphFrame = new SimpleInternalFrame("Workflow view", graphPanel.getButtonToolBar(), graphPanel.getScroller());
		}
		return graphFrame;
	}

	public SimpleInternalFrame getTreeFrame() throws MicroarrayException {
		if (treeFrame == null) {
			treeFrame = new SimpleInternalFrame("Datasets");
			treeFrame.add(tree);
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
			addPropertyChangeListener(detailsFrameWithListener);

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
			} catch (UnsupportedLookAndFeelException e1) {
				e1.printStackTrace();
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

		if (folderName == null) {
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

						String dataSetName = ImportUtils.convertToDatasetName(item.getOutput().getName());
						ContentType contentType = item.getType();
						Object dataSource = item.getInput();


						// Selects folder where data is imported to, or creates a
						// new one
						DataFolder folder = initializeFolderForImport(folderName);

						// get the InputStream for the data source
						InputStream input;

						if (dataSource instanceof File) {
							input = new FileInputStream((File) (dataSource));

						} else if (dataSource instanceof URL) {
							// TODO Not used anymore, URL-files are saved to the
							// temp file
							URL url = (URL) dataSource;
							try {
								input = url.openStream();

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

						// create new data
						DataBean data = manager.createDataBean(dataSetName, input);
						data.setContentType(contentType);

						// add the operation (all databeans have their own import
						// operation
						// instance, it would be nice if they would be grouped)
						Operation importOperation = new Operation(OperationDefinition.IMPORT_DEFINITION, new DataBean[] { data });
						data.setOperation(importOperation);

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
					// select data
					final DataBean selectedBean = lastGroupMember;
					SwingUtilities.invokeAndWait(new Runnable() {
						public void run() {
							getSelectionManager().selectSingle(selectedBean, this);
						}
					});

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
		showDialog(title, message, details, severity, modal, ChipsterDialog.DetailsVisibility.DETAILS_HIDDEN);
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
	public void showDialog(String title, String message, String details, Severity severity, boolean modal, ChipsterDialog.DetailsVisibility detailsVisibility) {
		DialogInfo dialogInfo = new DialogInfo(severity, title, message, details);
		ChipsterDialog.showDialog(mainFrame, dialogInfo, detailsVisibility, modal);
	}

	@Override
	public File saveWorkflow() {

		try {
			workflowManager.initialiseScriptDirectory();
			JFileChooser fileChooser = this.getWorkflowFileChooser();
			int ret = fileChooser.showSaveDialog(this.getMainFrame());
			if (ret == JFileChooser.APPROVE_OPTION) {
				File selected = fileChooser.getSelectedFile();
				File newFile = selected.getName().endsWith(WorkflowManager.SCRIPT_EXTENSION) ? selected : new File(selected.getCanonicalPath() + "." + WorkflowManager.SCRIPT_EXTENSION);

				workflowManager.saveSelectedWorkflow(newFile);
				unsavedChanges = false;
				menuBar.updateMenuStatus();
				return newFile;
			}
			menuBar.updateMenuStatus();

		} catch (IOException e) {
			reportException(e);
		}
		return null;
	}

	@Override
	public File openWorkflow() {

		try {
			workflowManager.initialiseScriptDirectory();
			JFileChooser fileChooser = this.getWorkflowFileChooser();
			int ret = fileChooser.showOpenDialog(this.getMainFrame());
			if (ret == JFileChooser.APPROVE_OPTION) {
				runWorkflow(fileChooser.getSelectedFile());
				unsavedChanges = false;
				menuBar.updateMenuStatus();
				return fileChooser.getSelectedFile();
			}
			menuBar.updateMenuStatus();

		} catch (IOException e) {
			reportException(e);
		}
		return null;
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
		String title = "Running " + task.getNamePrettyPrinted() + " failed. ";
		String message = "You may have used a tool or parameters which are unsuitable for the selected dataset, or\n" + "there might be a bug in the analysis tool itself.\n\n" + "The details below may provide hints about the problem. The most useful information is usually at the few last lines.";
		DialogInfo dialogInfo = new DialogInfo(Severity.WARNING, title, message, task.getScreenOutput());
		ChipsterDialog.showDialog(mainFrame, dialogInfo, ChipsterDialog.DetailsVisibility.DETAILS_ALWAYS_VISIBLE, false);
	}

	public void reportException(Exception e) {

		// collect error information to dialogInfo
		DialogInfo dialogInfo = new DialogInfo(Severity.ERROR, "An error has occurred and the action was not performed successfully.", "If problem persist, please check that your data is valid. For more information open the details panel below.", null);

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
		ChipsterDialog.showDialog(this.mainFrame, dialogInfo, ChipsterDialog.DetailsVisibility.DETAILS_HIDDEN, false);

		// we'll always output these to console and log for traceability and
		// easier IDE navigation
		e.printStackTrace();
		if (logger != null) {
			logger.error("client got exception", e);
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
				return VisualisationMethod.getDefaultVisualisationFor(beans.get(0));
			} else if (beans.size() > 1)
				for (VisualisationMethod method : VisualisationMethod.orderedDefaultCandidates()) {
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

			File file = ImportUtils.createTempFile(ImportUtils.convertToDatasetName(ImportUtils.URLToFilename(selectedURL)), ImportUtils.getExtension(ImportUtils.URLToFilename(selectedURL)));

			ImportUtils.getURLFileLoader().loadFileFromURL(selectedURL, file, importFolder, urlImportDlg.isSkipSelected());
		}
	}

	public void openClipboardImport() throws MicroarrayException, IOException {
		new ClipboardImportDialog(this);
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
					saveSession();
				} catch (Exception exp) {
					this.showErrorDialog("Session saving failed", exp);
					return;
				}

				// cancel also if the dialog is closed from the corner x
			} else if (returnValue != 1) {
				return;
			}
		}

		// hide immediately to look more reactive...
		childScreens.disposeAll();
		mainFrame.setVisible(false);

		super.quit();

		// this closes the application
		mainFrame.dispose();
		System.exit(0);
	}

	public static void start() {
		ClientListener shutdownListener = new ClientListener() {
			public void onSuccessfulInitialisation() {
				// do nothing
			}

			public void onFailedInitialisation() {
				System.exit(1);

			}
		};

		start(shutdownListener, null);
	}

	/**
	 * Starts Chipster client. Configuration (logging) should be initialised
	 * before calling this method.
	 * 
	 * @param shutdownListener
	 * 
	 */
	public static void start(final ClientListener clientListener, AuthenticationRequestListener overridingARL) {
		try {
			new SwingClientApplication(clientListener, overridingARL);
		} catch (Throwable t) {
			t.printStackTrace();
			if (logger != null) {
				logger.error(t.getMessage());
				logger.error(t);
			}
		}
	}

	public static void main(String[] args) throws IOException {
		try {
			MicroarrayConfiguration.loadConfiguration();

		} catch (OldConfigurationFormatException e) {
			reportOldConfigurationFormatException(e);
		}

		start();
	}

	public static void reportOldConfigurationFormatException(OldConfigurationFormatException e) {
		DialogInfo dialogInfo = new DialogInfo(Severity.ERROR, "Chipster configuration file not compatible", "Configuration file is from the previous version of the software. You can remove nami-work-files directory from your home directory so the new configuration is created automatically.", "Reason: " + e.getMessage());
		ChipsterDialog.showDialog(null, dialogInfo, DetailsVisibility.DETAILS_HIDDEN, true);
		throw new RuntimeException("configuration not compatible, will not start");
	}

	@Override
	public void showSourceFor(String operationName) throws TaskException {
		childScreens.show("ShowSource", true, operationName);
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
	protected AuthenticationRequestListener getAuthenticationRequestListener() {

		AuthenticationRequestListener authenticator;

		if (overridingARL != null) {
			authenticator = overridingARL;
		} else {
			authenticator = new Authenticator();
		}

		authenticator.setLoginListener(new ClientLoginListener() {
			public void firstLogin() {
				try {
					initialiseGUI();
				} catch (Exception e) {
					reportException(e);
				}
			}

			public void loginCancelled() {
				System.exit(1);
			}
		});

		return authenticator;
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

		// FIXME All files to default filter and other as a single file
		// filters
		String description = "Common microarray filetypes (cel, spot, gpr, txt, csv, tsv)";
		String[] extensions = { "cel", // affymetrix
		"spot", // SPOT files
		"gpr", // GenePix
		"txt", "csv", // illumina
		"tsv" // chipster
		};

		for (FileFilter filter : importExportFileChooser.getChoosableFileFilters()) {
			importExportFileChooser.removeChoosableFileFilter(filter);
		}

		importExportFileChooser.setFileSelectionMode(JFileChooser.FILES_ONLY);
		importExportFileChooser.setAcceptAllFileFilterUsed(true);
		FileFilter filter = new GeneralFileFilter(description, extensions);
		importExportFileChooser.addChoosableFileFilter(filter);
		importExportFileChooser.setFileFilter(filter);

		ImportSettingsAccessory access = new ImportSettingsAccessory(importExportFileChooser);
		importExportFileChooser.setAccessory(access);

		fixFileChooserFontSize(importExportFileChooser);

		return importExportFileChooser;
	}

	private JFileChooser getSnapshotFileChooser(JComponent accessory) {
		if (snapshotFileChooser == null) {
			snapshotFileChooser = ImportUtils.getFixedFileChooser();

			String[] extensions = { FSSnapshottingSession.SNAPSHOT_EXTENSION };
			snapshotFileChooser.setFileFilter(new GeneralFileFilter("Chipster Session", extensions));
			snapshotFileChooser.setSelectedFile(new File("session." + FSSnapshottingSession.SNAPSHOT_EXTENSION));
			snapshotFileChooser.setAcceptAllFileFilterUsed(false);
			snapshotFileChooser.setMultiSelectionEnabled(false);

		}
		snapshotFileChooser.setAccessory(accessory);

		fixFileChooserFontSize(snapshotFileChooser);

		return snapshotFileChooser;
	}

	private JFileChooser getWorkflowFileChooser() {
		if (workflowFileChooser == null) {
			workflowFileChooser = ImportUtils.getFixedFileChooser(WorkflowManager.SCRIPT_DIRECTORY);

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
	 * @param file
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
				}
			}
		});
		waitPanel.startWaiting("Please wait while " + taskName + "...");
		backgroundThread.start();
	}

	public void viewHelpFor(OperationDefinition definition) {
		viewHelp(HelpMapping.mapToHelppage(definition));
	}

	public void viewHelp(String page) {
		try {
			BrowserLauncher.openURL("https://extras.csc.fi/biosciences/" + page);

		} catch (Exception e) {
			reportException(e);
		}
	}

	protected void garbageCollect() {
		System.gc();
		statusBar.updateMemoryIndicator();
	}

	public TaskManagerScreen getTaskManagerScreen() {
		TaskExecutor jobExecutor = Session.getSession().getJobExecutor("client-job-executor");
		return new TaskManagerScreen(jobExecutor);
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
	public void loadSession() {

		SnapshotAccessory accessory = new SnapshotAccessory();
		final JFileChooser fileChooser = getSnapshotFileChooser(accessory);
		int ret = fileChooser.showOpenDialog(this.getMainFrame());
		final ClientApplication application = this; // for inner class
		
		if (ret == JFileChooser.APPROVE_OPTION) {
				if (accessory.clearSession()) {
					if (!clearSession()) {
						return; // loading cancelled
					}
				}
				
				runBlockingTask("loading session", new Runnable() {
					public void run() {						
						try {
							final List<DataItem> newItems = manager.loadSnapshot(fileChooser.getSelectedFile(), manager.getRootFolder(), application);
							SwingUtilities.invokeAndWait(new Runnable() {
								public void run() {
									getSelectionManager().selectSingle(newItems.get(newItems.size() - 1), this); // select last
								}
							});
						} catch (Exception e) {
							throw new RuntimeException(e);
						}						
					}
				});
		}
		menuBar.updateMenuStatus();
		unsavedChanges = false;
	}

	@Override
	public void saveSession() {

		JFileChooser fileChooser = getSnapshotFileChooser(null);
		int ret = fileChooser.showSaveDialog(this.getMainFrame());
		final ClientApplication application = this; // for inner class
		
		if (ret == JFileChooser.APPROVE_OPTION) {
			try {
				final File file = fileChooser.getSelectedFile().getName().endsWith("." + FSSnapshottingSession.SNAPSHOT_EXTENSION) ? fileChooser.getSelectedFile() : new File(fileChooser.getSelectedFile().getCanonicalPath() + "." + FSSnapshottingSession.SNAPSHOT_EXTENSION);

				if (file.exists()) {
					int returnValue = JOptionPane.DEFAULT_OPTION;

					String message = "The file " + file.getCanonicalPath() + " exists already. Do you want " + "to replace it with the one you are saving?";

					Object[] options = { "Cancel", "Replace" };

					returnValue = JOptionPane.showOptionDialog(this.getMainFrame(), message, "Confirm replace", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE, null, options, options[0]);

					if (returnValue != 1) {
						return;
					}

				}

				runBlockingTask("saving session", new Runnable() {

					public void run() {

						try {
							getDataManager().saveSnapshot(file, application);
						} catch (IOException e) {
							throw new RuntimeException(e);
						}
					}
				});
				
				menuBar.updateMenuStatus();
				setUnsavedChanges(false);
				
			} catch (Exception exp) {
				showErrorDialog("Saving session failed.", exp);
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
	public void heartBeat() {
		statusBar.updateMemoryIndicator();
	}

	public DataManager getDataManager() {
		return manager;
	}

	public void setUnsavedChanges(boolean value) {
		this.unsavedChanges = value;
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
}