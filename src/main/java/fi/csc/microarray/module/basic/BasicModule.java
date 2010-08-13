package fi.csc.microarray.module.basic;

import java.util.List;

import javax.swing.JMenu;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.TypeTag;
import fi.csc.microarray.databeans.features.RestrictModifier;
import fi.csc.microarray.databeans.features.stat.LogModifier;
import fi.csc.microarray.databeans.features.stat.NegModifier;
import fi.csc.microarray.databeans.features.table.RowCountProvider;
import fi.csc.microarray.databeans.features.table.TableColumnProvider;
import fi.csc.microarray.module.Module;

public class BasicModule implements Module {

	public void plugContentTypes(DataManager manager) {

		manager.plugContentType("text/plain", false, false, "plain text", VisualConstants.ICON_TYPE_TEXT, "txt", "dat", "wee");
		manager.plugContentType("application/octet-stream", false, true, "binary", VisualConstants.ICON_TYPE_BINARY, "");
		
		manager.plugContentType("text/tab", false, false, "tab separated values", VisualConstants.ICON_TYPE_TABLE, "tsv");
		manager.plugContentType("text/csv", false, false, "comma separated values", VisualConstants.ICON_TYPE_TABLE, "csv");

		manager.plugContentType("image/png", true, true, "PNG image", VisualConstants.ICON_TYPE_IMAGE, "png");
		manager.plugContentType("image/gif", true, true, "GIF image", VisualConstants.ICON_TYPE_IMAGE, "gif");
		manager.plugContentType("image/jpeg", true, true, "JPEG image", VisualConstants.ICON_TYPE_IMAGE, "jpg", "jpeg");
		
		manager.plugContentType("text/html", true, false, "HTML document", VisualConstants.ICON_TYPE_HTML, "html", "htm");
	}

	public void plugFeatures(DataManager manager) {
		manager.plugFeatureFactory("/column", new TableColumnProvider());
		manager.plugFeatureFactory("/rowcount", new RowCountProvider());
	}

	public void plugModifiers(DataManager manager) {
		manager.plugModifier("log", new LogModifier());
		manager.plugModifier("neg", new NegModifier());
		manager.plugModifier("restrict", new RestrictModifier());
	}

	@Override
	public void plugTypeTags(DataManager manager) {
		manager.plugTypeTag(new TypeTag("table-without-header"));
		
	}

	@Override
	public String getServerModuleName() {
		return null;
	}

	@Override
	public void addImportMenuItems(JMenu importMenu) {
		// do nothing
	}

	@Override
	public void addImportLinks(QuickLinkPanel quickLinkPanel, List<JXHyperlink> importLinks) {
		// do nothing
	}
}
