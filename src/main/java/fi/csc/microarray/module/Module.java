package fi.csc.microarray.module;

import java.util.List;

import javax.swing.JMenu;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.databeans.DataManager;

public interface Module {

	public void plugFeatures(DataManager manager);
	public void plugModifiers(DataManager manager);
	public void plugContentTypes(DataManager manager);
	public void plugTypeTags(DataManager manager);
	public String getServerModuleName();
	public void addImportMenuItems(JMenu importMenu);
	public void addImportLinks(QuickLinkPanel quickLinkPanel, List<JXHyperlink> importLinks);
}
