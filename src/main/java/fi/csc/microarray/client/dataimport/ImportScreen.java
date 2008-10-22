package fi.csc.microarray.client.dataimport;

import java.awt.BorderLayout;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.io.File;
import java.io.IOException;
import java.nio.channels.ClosedByInterruptException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSplitPane;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.SwingClientApplication;
import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.dataimport.table.InformationDialog;
import fi.csc.microarray.client.dataimport.table.TableInternalFrame;
import fi.csc.microarray.client.dataimport.tools.ToolsInternalFrame;
import fi.csc.microarray.client.dataimport.trimmer.DataTrimmer;
import fi.csc.microarray.client.screen.ScreenBase;

/**
 * The Import tool main class. Includes the GUI and also functions which 
 * gives access to the classes which do the file reading/writing etc.
 * 
 * @author mkoski klemela
 *
 */
public class ImportScreen extends ScreenBase  implements ImportScreenModel, ActionListener, WindowListener{
	
	private static final int TOOLS_FRAME_WIDTH = 200;
	
	private static final Dimension BUTTON_SIZE = new Dimension(80, 22);
	private static final Dimension IMPORT_SCREEN_SIZE = new Dimension(850,700);
	
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(ImportScreen.class);
		
	private JFrame frame;
	private JSplitPane mainSplit;
	private ConversionModel conversionModel;
	private ToolsInternalFrame toolsFrame;
	private TableInternalFrame tableFrame;
	private JPanel changeStepPanel;
	
	/**
	 * Enumeration for import steps
	 */
	public enum Step {
		FIRST, 
		SECOND 
	};
	
	/**
	 * Current step
	 */
	private Step currentStep;
	
	/**
	 * Class which takes care of column types and chip number on the second step
	 */
	private ColumnTypeManager columnTypeManager;
	
	private JButton helpButton;
	private JButton cancelButton;			
	private JButton backButton;
	private JButton nextButton;
	private JButton finishButton;
	
	private ClientApplication application;
	private DataTrimmer dataTrimmer;
	private DataTrimmer flagTrimmer;

	private ImportSession importSession;
	
	private Iterator<File> files;
	
	/**
	 * Helper class to write files after correct file description are set
	 * 
	 * @author mkoski
	 *
	 */
	class WriteToFileProcess extends RunnableImportProcess{

		public WriteToFileProcess(ProgressInformator informator) {
			super(informator);
		}

		@Override
		public void taskToDo() {
			List<ImportItem> dataItems = new ArrayList<ImportItem>();
			
			boolean hasNext = files.hasNext();
			File currentFile = null;
			try {
				do{

					// Write file to disk
					conversionModel.writeToFile(informator);		// Throws exceptions

					// Import files to application and create the DataBeans
					if (application != null) {
						
						
						for(File outputFile : conversionModel.getOutputFiles()){
							ImportItem item = new ImportItem(outputFile);
							item.setType(Session.getSession().getDataManager().getContentType("text/tab"));
							item.setFilename(outputFile.getName());
							
							dataItems.add(item);
						}												
					}

					// Get next file
					hasNext = files.hasNext();
					if(hasNext){
						currentFile = files.next();
						conversionModel.setInputFile(currentFile);
					}

					if(hasNext){
						informator.setMessage("Writing data to disk from file " + currentFile.getName());
					}

				} while(hasNext && importSession.getUseSameDescriptions());
				
				application.importGroup(dataItems, importSession.getDestinationFolder());
			}
			
			// Catch exceptions
			catch (ClosedByInterruptException cbie){
				// Do nothing
			}
			catch (Exception ex) {
				if (application != null) {
					application.reportException(ex);
				} else {
					ex.printStackTrace();
				}
			} finally{
				
				// All files wrote with the same description
				if(!hasNext){
					// Reset the import screen and do not show it after resetting
					resetImportScreen(false);
					resetImportSession();
					
					// Initialize the first step, but do no show it
					gotoFirstStep();
				} else {
					// Reset the import screen and show it again after resetting
					resetImportScreen(true);
					
					conversionModel.setInputFile(currentFile);
					
					// Update table and chop the data from input file
					updateTable(false);
					
					// Initialize the first step and show it
					gotoFirstStep();
				}
			}
		}
	}
	
	/**
	 * Helper class to run import process
	 * 
	 * @author mkoski
	 *
	 */
	class UpdateTableProcess extends RunnableImportProcess{
		
		/**
		 * Parse headers and footers? On the first step this should be true 
		 * and on the second step false
		 */
		private boolean hideHeaderFooter;
		
		private Object[][] choppedDataMatrix;
		private String[] columnTitles;
		
		
		public UpdateTableProcess(ProgressInformator informator, boolean hideHeaderFooter) {
			super(informator);
			this.hideHeaderFooter = hideHeaderFooter;
		}

		public void taskToDo() {
			logger.debug("Starting import process");
			
			choppedDataMatrix = null;
			try {
				choppedDataMatrix = conversionModel.chopData(hideHeaderFooter, this.getInformator());
			} catch (IOException e) {
				application.reportException(e);
			}
			
			columnTitles = conversionModel.getColumnTitles();
			
			tableFrame.getTable().setData(UpdateTableProcess.this.choppedDataMatrix, UpdateTableProcess.this.columnTitles);			
			
			
			SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					// To show right column delimeter selected
					ImportScreen.this.getToolsFrame().updateDelimeterPanel();	
					
					// Updates "Showing columns..." label
					ImportScreen.this.getTableFrame().updateShowingColumnsLabel();
				}
			} );
			
			
			logger.debug("Data set to table");
		}
	}
	
	/**
	 * Constructor method. This method uses resetImportScreen to create all 
	 * necessery components. It also sets look'n'feel settings and 
	 * gets ready to show the first import step
	 *
	 */
	public ImportScreen(){
		// Reset import screen, but do not show it after resetting
		resetImportScreen(false);
		resetImportSession();
		SwingClientApplication.setPlastic3DLookAndFeel(frame);
		//updateUI();
		gotoFirstStep();
		this.getFrame().addWindowListener(this);
	}

	/**
	 * Updates the internal frames UI after look'n'feel is changed.
	 *
	 *
	private void updateUI() {
		toolsFrame.updateUI();
		tableFrame.updateUI();
	}*/

	/**
	 * Panel which includes cancel, back, next and finish buttons
	 * 
	 * @return panel for cancel, back, next and finish buttons
	 */
	private Component getChangeStepButtonsPanel() {
		if( changeStepPanel == null){
			changeStepPanel = new JPanel(new GridBagLayout());
			
			//TODO icon
			helpButton = new JButton("Help");
			cancelButton = new JButton("Cancel", VisualConstants.IMPORT_CANCEL_ICON);			
			backButton = new JButton("Back", VisualConstants.IMPORT_BACK_ICON);
			nextButton = new JButton("Next", VisualConstants.IMPORT_NEXT_ICON);
			finishButton = new JButton("Finish", VisualConstants.IMPORT_FINISH_ICON);
			nextButton.setHorizontalTextPosition(SwingConstants.LEADING);
			finishButton.setHorizontalTextPosition(SwingConstants.LEADING);
			cancelButton.setHorizontalTextPosition(SwingConstants.LEADING);
			
			helpButton.setPreferredSize(BUTTON_SIZE);
			cancelButton.setPreferredSize(BUTTON_SIZE);
			backButton.setPreferredSize(BUTTON_SIZE);
			nextButton.setPreferredSize(BUTTON_SIZE);
			finishButton.setPreferredSize(BUTTON_SIZE);			
			
			helpButton.addActionListener(this);
			cancelButton.addActionListener(this);
			backButton.addActionListener(this);
			nextButton.addActionListener(this);
			finishButton.addActionListener(this);
			
			GridBagConstraints c = new GridBagConstraints();
			c.weightx = 0;
			c.fill= GridBagConstraints.NONE;
			c.insets.set(10,10,10,0);
			changeStepPanel.add(helpButton,c);
			c.weightx = 1;
			c.fill = GridBagConstraints.HORIZONTAL;			
			changeStepPanel.add(new JLabel(),c);
			c.weightx = 0;
			c.fill= GridBagConstraints.NONE;
			c.insets.set(10,10,10,0);
			changeStepPanel.add(backButton,c);
			changeStepPanel.add(nextButton,c);
			changeStepPanel.add(finishButton,c);
			c.insets.set(10,30,10,10);
			changeStepPanel.add(cancelButton,c);
		}
		return changeStepPanel;
	}
	
	/**
	 * Gets current step
	 * @return
	 */
	public Step getCurrentStep(){
		return currentStep;
	}
	
	/**
	 * Goes to the first import step. Calls methods which initialized 
	 * the table and other components for the first step.
	 */
	private void gotoFirstStep(){
		//This should be first to get updates right
		currentStep = Step.FIRST;
		
		// Removes the internal frames while updating
		this.getFrame().remove(tableFrame);
		this.getFrame().remove(toolsFrame);
		
		// Order is siqnificant in avoiding nullPointers
		tableFrame.initializeFirstStep();
		toolsFrame.initializeFirstStep();
		
		backButton.setEnabled(false);
		nextButton.setEnabled(true);
		finishButton.setEnabled(false);
	}
	
	/**
	 * Goes to the second import step. Calls methods which initialized 
	 * the table and other components for the second step.
	 */
	private void gotoSecondStep(){
		//This should be first to get updates right
		currentStep = Step.SECOND;
		
		// Order is siqnificant in avoiding nullPointers
		tableFrame.initializeSecondStep();
		toolsFrame.initializeSecondStep();
		
		backButton.setEnabled(true);		
		nextButton.setEnabled(false);
		finishButton.setEnabled(true);
		
	}
	
	public void actionPerformed(ActionEvent e){
		Object source = e.getSource();
		
		if (source == helpButton){
			if (this.getCurrentStep() == Step.FIRST){
				application.viewHelp("chipster-manual/import-help.html#step1");
			} else {
				application.viewHelp("chipster-manual/import-help.html#step2");
			}
		}
		
		// Back to first step
		if (source == backButton) {
			this.gotoFirstStep();
			this.updateTable(false);
		} 
		
		// To the second step
		else if (source == nextButton) {
			this.gotoSecondStep();
			this.updateTable(true);
		} 
		
		// Finish, write the data
		else if (source == finishButton) {
			this.finishButtonPressed();
			
		} 
		
		else if (source == cancelButton) {
			resetImportScreen(false);
			resetImportSession();
			gotoFirstStep();
			frame.setVisible(false);
		}
	}

	public boolean hasFrame() {
		return frame != null;
	}
	
	public JFrame getFrame() {
		return frame;
	}
	
	public ConversionModel getConversionModel(){
		return this.conversionModel;
	}
	
	public void setImportSession(ImportSession importSession){
		this.importSession = importSession;
		this.files = importSession.getCustomFiles().iterator();
		
		// Set the first input file
		conversionModel.setInputFile(files.next());
	}
	
	/**
	 * Updates the table and does file reading and chopping at the same time. 
	 * This operation is time consuming and that's why it is run on its own 
	 * thread provided by <code>RunnableImportProcess</code>
	 * 
	 */
	public void updateTable(boolean hideHeaderFooter) {
		ProgressInformator infoDialog = new InformationDialog("Updating table", "", this);
		UpdateTableProcess process = new UpdateTableProcess(infoDialog, hideHeaderFooter);
		process.runProcess();
	}

	public TableInternalFrame getTableFrame() {
		return tableFrame;
	}

	public ToolsInternalFrame getToolsFrame() {
		return toolsFrame;
	}

	public ColumnTypeManager getColumnTypeManager() {
		return columnTypeManager;
	}

	public DataTrimmer getDataTrimmer() {
		return this.dataTrimmer;
	}

	public DataTrimmer getFlagTrimmer() {
		return this.flagTrimmer;
	}
	
	/**
	 * Reset import screen. This can be used while creating totally new instance of 
	 * import screen or just resetting the default.
	 * 
	 * @param setVisibleAfterReset set import screen frame visible after 
	 *        it has been reseted
	 */
	public void resetImportScreen(boolean setVisibleAfterReset){
		
		application = Session.getSession().getApplication();
		
		// If frame exists, set visibility to false while resetting
		if(frame != null){
			frame.setVisible(false);
		}
		
		// Create new instances of conversion model, column type manager and trimmers
		conversionModel = new ConversionModel(this);
		columnTypeManager = new ColumnTypeManager(0);
		dataTrimmer = new DataTrimmer();
		
		// Ignore first
		dataTrimmer.addIgnoreColumnNumber(0);
		
		flagTrimmer = new DataTrimmer();
		flagTrimmer.addIgnoreColumnNumber(0);
		
		// Table has to be first to get references right
		tableFrame = new TableInternalFrame(this);
		toolsFrame = new ToolsInternalFrame(this);
		
		conversionModel.addConversionChangeListener(tableFrame);
		columnTypeManager.addColumnTypeChangeListener(tableFrame);
		columnTypeManager.addColumnTypeChangeListener(toolsFrame);
		
		frame = new JFrame("Import tool");
		mainSplit = new JSplitPane();
		
		frame.setLayout(new BorderLayout());
		frame.setSize(IMPORT_SCREEN_SIZE);
		
		mainSplit.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
		mainSplit.setDividerLocation(TOOLS_FRAME_WIDTH);
		mainSplit.setLeftComponent(toolsFrame);
		mainSplit.setRightComponent(tableFrame);
		mainSplit.getLeftComponent().setMinimumSize(new Dimension(150, 0));
		mainSplit.setResizeWeight(0);
		
		frame.add(mainSplit, BorderLayout.CENTER);
		
		// Reset buttons
		changeStepPanel = null;
		frame.add(getChangeStepButtonsPanel(), BorderLayout.SOUTH);
		
		if(setVisibleAfterReset){
			frame.setVisible(true);
		}
	}
	
	private void finishButtonPressed(){
		// No chips selected
		if(columnTypeManager.getChipCount() == 0){
			JOptionPane.showMessageDialog(this.getFrame(), "No chips to import. Select at least one sample.", "No chips selected", JOptionPane.ERROR_MESSAGE);
		} 
		
		// Write data to file
		else {
			RunnableImportProcess writeToFileProcess = new WriteToFileProcess(new InformationDialog("Writing data to disk", "Writing data to disk from file " + conversionModel.getInputFileName(), this));
			writeToFileProcess.runProcess();
		}
	}

	/**
	 * Resets import screen
	 */
	public void windowClosing(WindowEvent e) {
		resetImportScreen(false);
		resetImportSession();
		gotoFirstStep();
	}
	
	private void resetImportSession() {
		importSession = null;
	}

	// Other WindowListener methods. No need to implement
	public void windowActivated(WindowEvent e) { }
	public void windowClosed(WindowEvent e) { }
	public void windowDeactivated(WindowEvent e) { }
	public void windowDeiconified(WindowEvent e) { }
	public void windowIconified(WindowEvent e) { }
	public void windowOpened(WindowEvent e) { }
}

