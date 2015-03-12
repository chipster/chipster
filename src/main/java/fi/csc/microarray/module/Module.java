package fi.csc.microarray.module;

import java.io.IOException;
import java.net.MalformedURLException;
import java.util.List;

import javax.jms.JMSException;
import javax.swing.Icon;
import javax.swing.JMenu;
import javax.swing.JPanel;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.DbSession;

/**
 * Client side module. Encapsulates all application area specific logic, e.g., DNA microarray module
 * encapsulates all DNA microarray specific functionality.
 * 
 * @author Aleksi Kallio
 *
 */
public interface Module {

	/**
	 * Plugs features of this module to given data manager.
	 * 
	 * @param manager data manager to plug into
	 */
	public void plugFeatures(DataManager manager);
	
	/**
	 * Plugs modifiers of this module to given data manager.
	 * 
	 * @param manager data manager to plug into
	 */
	public void plugModifiers(DataManager manager);
	
	/**
	 * Plugs content types of this module to given data manager.
	 * 
	 * @param manager data manager to plug into
	 */
	public void plugContentTypes(DataManager manager);
	
	/**
	 * Returns the name of the server module associated to this 
	 * this client side module, or null if not available.
	 *  
	 * @return String array of server module names, or null
	 */
	public String[] getServerModuleNames();
	
	/**
	 * Adds import menu items to menu.
	 * 
	 * @param importMenu menu to add links to
	 */
	public void addImportMenuItems(JMenu importMenu);
	
	/**
	 * Adds import links link list.
	 * 
	 * @param quickLinkPanel quick link panel used when creating the links
	 * @param importLinks link list to add to
	 */
	public void addImportLinks(QuickLinkPanel quickLinkPanel, List<JXHyperlink> importLinks);
	
	/**
	 * Is import tool supported when this module is used? 
	 * 
	 * @return true iff import tool is supported by this module
	 */
	public boolean isImportToolSupported();
	
	/**
	 * Checks data bean for workflow compatibility.
	 * 
	 * @param data data bean to check
	 * @return true iff data bean supports workflows under, when this module is used
	 */
	public boolean isWorkflowCompatible(DataBean data);
	
	/**
	 * Checks if data bean is metadata. Metadata
	 * datasets are not considered to be part of the normal 
	 * workflow ancestor relations, but instead they 
	 * are considered to be annotations that are linked
	 * to real datasets.
	 * 
	 * @param data data bean to check
	 * @return true iff data bean is metadata 
	 */
	public boolean isMetadata(DataBean data);
	
	/**
	 * Pre-processes metadata datasets before the task is started.
	 * 
	 * @param operation operation that is going to run 
	 * @param metadataOutput the metadata input to pre-process
	 * @throws MicroarrayException
	 * @throws IOException
	 */
	void preProcessInputMetadata(Operation oper, DataBean metadataInput) throws MicroarrayException, IOException;
	
	/**
	 * Post-processes metadata datasets after they have been created by finished tasks.
	 * 
	 * @param operationRecord operation that generated the output metadata 
	 * @param metadataOutput the metadata output to post-process
	 * @throws MicroarrayException
	 * @throws IOException
	 */
	public void postProcessOutputMetadata(OperationRecord operationRecord, DataBean metadataOutput) throws MicroarrayException, IOException;
	
	/**
	 * Returns visualisation methods of this module.
	 * 
	 * @return array of visualisation methods
	 */
	public VisualisationMethod[] getVisualisationMethods();
	
	/**
	 * Return URL to example session file or null, if not available.
	 * 
	 * @param isStandalone true iff client is running in standalone mode
	 * @return url or null
	 * @throws MalformedURLException
	 * @throws Exception 
	 * @throws JMSException 
	 */
	public List<DbSession> getExampleSessions(boolean isStandalone) throws JMSException, Exception;

	/**
	 * If module is bundled with a repository of workflows, returns them.
	 * Workflows are represented by a String array of two cells:
	 * display name and resource name of the bundled workflow file.
	 * 
	 * @return array of String arrays, zero length if none available
	 */
	public String[][] getRepositoryWorkflows();

	/**
	 * Read the name of the operation and return a short version suitable for 
	 * small UI components.
	 * 
	 * @param operation the category we want to shorten
	 * 
	 * @return short name (4 letters or less)
	 */
	public String getShortDataName(String categoryName);

	/**
	 * Should workflow engine check for the number of results? If the module contains tools
	 * that produce variable amounts of results, then number of results should not be checked.
	 * 
	 * @return should workflow engine check for the number of results?
	 */
	public boolean countOperationResults();
	
	public boolean notesVisibleAtStartup();

	/**
	 * Create context link panel that contains links for selecting, visualising etc. for datasets. 
	 * 
	 * @param selectedDataCount number of currently selected datasets 
	 * 
	 * @return panel or null if context link panel should not be shown
	 */
	public JPanel getContextLinkPanel(int selectedDataCount);

	/**
	 * Name of the module for the UI
	 * 
	 * @return
	 */
	public String getDisplayName();
	
	/**
	 * Url to manual home page.
	 * 
	 * @return
	 */
	public String getManualHome();

	/**
	 * Flags spreadsheet columns that support linking by this module.
	 *  
	 * @param columns spreadsheet columns
	 * @param data DataBean TagTypes are used for finding the right columns
	 * 
	 * @return Boolean list with the same size as columnNames
	 * @throws MicroarrayException 
	 */
	public List<Boolean> flagLinkableColumns(Table columns, DataBean data);

	/**
	 * Converts spreadshoot cell into linkable {@link IntegratedEntity}.
	 * 
	 * @param columns current row 
	 * @param column index of cell 
	 * @param data 
	 * 
	 * @return {@link IntegratedEntity} that is point selected when link is clicked
	 */
	public IntegratedEntity createLinkableEntity(Table columns, DataBean data);

	/**
	 * Converts server module name into GUI friendly name.
	 */
	public String getModuleLongName(String moduleName);

	/**
	 * Looks into data bean, possibly reading a little bit of content data, and adds type tags appropriate for the module.
	 * 
	 * @param data data bean to add type tags to
	 * @throws MicroarrayException 
	 * @throws IOException 
	 */
	public void addTypeTags(DataBean data) throws MicroarrayException, IOException;

	
	/**
	 * Returns icon for given data, depending on it's type and possibly content.
	 * 
	 * @return icon
	 */
	public Icon getIconFor(DataBean data);
	
}
