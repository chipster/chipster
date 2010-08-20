package fi.csc.microarray.module;

import java.util.List;

import javax.swing.JMenu;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;

/**
 * Client side module. Encapsulates all application area specific logic, e.g., DNA microarray module
 * encapsulates all DNA microarray specific functionality.
 * 
 * @author Aleksi Kallio
 *
 */
public interface Module {

	public void plugFeatures(DataManager manager);
	public void plugModifiers(DataManager manager);
	public void plugContentTypes(DataManager manager);
	public void plugTypeTags(DataManager manager);
	public String getServerModuleName();
	public void addImportMenuItems(JMenu importMenu);
	public void addImportLinks(QuickLinkPanel quickLinkPanel, List<JXHyperlink> importLinks);
	public boolean isImportToolSupported();
	public boolean isWorkflowCompatible(DataBean data);
	public VisualisationMethod[] getVisualisationMethods();
}
