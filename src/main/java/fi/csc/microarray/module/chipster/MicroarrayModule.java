package fi.csc.microarray.module.chipster;

import java.awt.event.ActionEvent;
import java.util.List;

import javax.swing.AbstractAction;
import javax.swing.JMenu;
import javax.swing.JMenuItem;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.TaskImportDialog;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.features.bio.EmbeddedBinaryProvider;
import fi.csc.microarray.databeans.features.bio.IdentifierProvider;
import fi.csc.microarray.databeans.features.bio.NormalisedExpressionProvider;
import fi.csc.microarray.databeans.features.bio.PhenodataProvider;
import fi.csc.microarray.databeans.features.stat.HierarchicalClusterProvider;
import fi.csc.microarray.databeans.features.stat.SomClusterProvider;
import fi.csc.microarray.module.Module;

public class MicroarrayModule implements Module {

	public static final String SERVER_MODULE_MICROARRAY = "microarray";

	public static final String ANNOTATION_ID = "annotate-genelist2html.R";

	public static final String IMPORT_FROM_ARRAYEXPRESS_ID = "import-from-ArrayExpress.R";
	public static final String IMPORT_FROM_GEO_ID = "import-soft2.R";

	public void plugContentTypes(DataManager manager) {
		manager.plugContentType("application/x-treeview", true, false, "Newick formatted tree from clustering", VisualConstants.ICON_TYPE_TEXT, "tre");
		manager.plugContentType("application/cel", true, false, "Affymetrix CEL", VisualConstants.ICON_TYPE_RAWDATA, "cel");
	}

	public void plugFeatures(DataManager manager) {
		manager.plugFeatureFactory("/normalised-expression", new NormalisedExpressionProvider());
		manager.plugFeatureFactory("/phenodata", new PhenodataProvider());
		manager.plugFeatureFactory("/identifier", new IdentifierProvider());
		manager.plugFeatureFactory("/embedded-binary-content", new EmbeddedBinaryProvider());
		manager.plugFeatureFactory("/clusters/som", new SomClusterProvider());
		manager.plugFeatureFactory("/clusters/hierarchical", new HierarchicalClusterProvider());
	}

	public void plugModifiers(DataManager manager) {
		// nothing to plug
	}

	@Override
	public void plugTypeTags(DataManager manager) {
		// TODO Auto-generated method stub

	}

	@Override
	public String getServerModuleName() {
		return SERVER_MODULE_MICROARRAY;
	}

	@Override
	public void addImportMenuItems(JMenu importMenu) {
		importMenu.add(getImportFromArrayExpressMenuItem());
		importMenu.add(getImportFromGEOMenuItem());
		importMenu.addSeparator();
	}

	private JMenuItem getImportFromArrayExpressMenuItem() {
		JMenuItem importFromArrayExpressMenuItem = new JMenuItem();
		importFromArrayExpressMenuItem.setText("ArrayExpress...");
		importFromArrayExpressMenuItem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				doImportFromArrayExpress();
			}
		});

		return importFromArrayExpressMenuItem;
	}

	private JMenuItem getImportFromGEOMenuItem() {
		JMenuItem importFromGEOMenuItem = new JMenuItem();
		importFromGEOMenuItem.setText("GEO...");
		importFromGEOMenuItem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				doImportFromGEO();
			}
		});
		return importFromGEOMenuItem;
	}

	@Override
	public void addImportLinks(QuickLinkPanel quickLinkPanel, List<JXHyperlink> importLinks) {

		importLinks.add(quickLinkPanel.createLink("Import from ArrayExpress ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				doImportFromArrayExpress();
			}
		}));

		importLinks.add(quickLinkPanel.createLink("Import from GEO ", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				doImportFromGEO();
			}
		}));
	}

	private void doImportFromGEO() {
		try {
			ClientApplication application = Session.getSession().getApplication();
			Operation importOperation = new Operation(application.getOperationDefinition(MicroarrayModule.IMPORT_FROM_GEO_ID), new DataBean[] {});
			new TaskImportDialog(application, "GEO", importOperation);
			
		} catch (Exception me) {
			Session.getSession().getApplication().reportException(me);
		}
	}

	private void doImportFromArrayExpress() {
		try {
			ClientApplication application = Session.getSession().getApplication();
			Operation importOperation = new Operation(application.getOperationDefinition(MicroarrayModule.IMPORT_FROM_ARRAYEXPRESS_ID), new DataBean[] {});
			new TaskImportDialog(application, "ArrayExpress", importOperation);
			
		} catch (Exception me) {
			Session.getSession().getApplication().reportException(me);
		}
	}

}
