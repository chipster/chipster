package fi.csc.microarray.module.chipster;

import java.awt.event.ActionEvent;
import java.net.MalformedURLException;
import java.net.URL;
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
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.methods.ArrayLayout;
import fi.csc.microarray.client.visualisation.methods.ClusteredProfiles;
import fi.csc.microarray.client.visualisation.methods.ExpressionProfile;
import fi.csc.microarray.client.visualisation.methods.HierarchicalClustering;
import fi.csc.microarray.client.visualisation.methods.Histogram;
import fi.csc.microarray.client.visualisation.methods.SOM;
import fi.csc.microarray.client.visualisation.methods.Scatterplot;
import fi.csc.microarray.client.visualisation.methods.Scatterplot3DPCA;
import fi.csc.microarray.client.visualisation.methods.VennDiagram;
import fi.csc.microarray.client.visualisation.methods.Volcanoplot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowser;
import fi.csc.microarray.client.visualisation.methods.threed.Scatterplot3D;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.TypeTag;
import fi.csc.microarray.databeans.features.bio.EmbeddedBinaryProvider;
import fi.csc.microarray.databeans.features.bio.IdentifierProvider;
import fi.csc.microarray.databeans.features.bio.NormalisedExpressionProvider;
import fi.csc.microarray.databeans.features.stat.HierarchicalClusterProvider;
import fi.csc.microarray.databeans.features.stat.SomClusterProvider;
import fi.csc.microarray.module.Module;

public class MicroarrayModule implements Module {

	private static final String EXAMPLE_SESSION_URL = "http://chipster.csc.fi/examples/ExampleSessionChipsterV2.cs";
	
	public static class TypeTags {
		public static final TypeTag RAW_AFFYMETRIX_EXPRESSION_VALUES  = new TypeTag(SERVER_MODULE_NAME, "raw-arrymetrix-expression-values", "must be in CEL format");;
		public static final TypeTag RAW_EXPRESSION_VALUES  = new TypeTag(SERVER_MODULE_NAME, "raw-expression-values", "");;
		public static final TypeTag NORMALISED_EXPRESSION_VALUES = new TypeTag(SERVER_MODULE_NAME, "normalised-expression-values", "must have columns following name pattern \"chip.*\"");
		public static final TypeTag SIGNIFICANT_EXPRESSION_FOLD_CHANGES = new TypeTag(SERVER_MODULE_NAME, "significant-expression-fold-changes", "must have columns following name patterns \"FC.*\" and \"p.*\"");
		public static final TypeTag EXPRESSION_PRIMARY_COMPONENTS_GENEWISE = new TypeTag(SERVER_MODULE_NAME, "expression-primary-components-genewise", "");
		public static final TypeTag EXPRESSION_PRIMARY_COMPONENTS_CHIPWISE = new TypeTag(SERVER_MODULE_NAME, "expression-primary-components-chipwise", "");
		public static final TypeTag ORDERED_GENOMIC_ENTITIES = new TypeTag(SERVER_MODULE_NAME, "ordered-genomic-entities", "");
		public static final TypeTag CLUSTERED_EXPRESSION_VALUES = new TypeTag(SERVER_MODULE_NAME, "clustered-expression-values", "must have column \"cluster\"");
		public static final TypeTag SOM_CLUSTERED_EXPRESSION_VALUES = new TypeTag(SERVER_MODULE_NAME, "som-clustered-expression-values", "must have columns \"colours\", \"distance2first\", \"cluster\", \"griddim\"");
	}
	
	public static class VisualisationMethods {
		public static VisualisationMethod ARRAY_LAYOUT = new VisualisationMethod("Array layout", ArrayLayout.class, VisualConstants.ARRAY_MENUICON, -1, 0.0009);
		public static VisualisationMethod HISTOGRAM = new VisualisationMethod("Histogram", Histogram.class, VisualConstants.HISTOGRAM_MENUICON, -1, 0.024);
		public static VisualisationMethod SCATTERPLOT = new VisualisationMethod("Scatterplot", Scatterplot.class, VisualConstants.SCATTER_MENUICON, -1, 0.039);
		public static VisualisationMethod SCATTERPLOT3D = new VisualisationMethod("3D Scatterplot", Scatterplot3D.class, VisualConstants.SCATTER3D_MENUICON, -1, 0.082);
		public static VisualisationMethod SCATTERPLOT3DPCA = new VisualisationMethod("3D Scatterplot for PCA", Scatterplot3DPCA.class, VisualConstants.SCATTER3DPCA_MENUICON, -1, 0.082);
		public static VisualisationMethod VOLCANOPLOT = new VisualisationMethod("Volcano plot", Volcanoplot.class, VisualConstants.VOLCANO_MENUICON, -1, 0.039); 
		public static VisualisationMethod SOM = new VisualisationMethod("SOM", SOM.class, VisualConstants.SOM_MENUICON, 3, 0.034);
		public static VisualisationMethod HIERARCHICAL = new VisualisationMethod("Hierarchical clustering", HierarchicalClustering.class, VisualConstants.HC_MENUICON, 3, 0.09);
		public static VisualisationMethod EXPRESSION_PROFILE = new VisualisationMethod("Expression profile", ExpressionProfile.class, VisualConstants.PROFILE_MENUICON, -1, 0.1);
		public static VisualisationMethod CLUSTERED_PROFILES = new VisualisationMethod("Clustered profiles", ClusteredProfiles.class, VisualConstants.PROFILES_MENUICON, -1, 0.087);
		public static VisualisationMethod VENN_DIAGRAM = new VisualisationMethod("Venn-diagram", VennDiagram.class, VisualConstants.VENN_MENUICON, 1, -1);
		public static VisualisationMethod GBROWSER = new VisualisationMethod("Genome browser", GenomeBrowser.class, VisualConstants.SCATTER_MENUICON, 1, -1);
	}
	
	public static final String SERVER_MODULE_NAME = "microarray";

	public static final String ANNOTATION_ID = "annotate-genelist2html.R";

	public static final String IMPORT_FROM_ARRAYEXPRESS_ID = "import-from-ArrayExpress.R";
	public static final String IMPORT_FROM_GEO_ID = "import-soft2.R";

	public void plugContentTypes(DataManager manager) {
		manager.plugContentType("application/x-treeview", true, false, "Newick formatted tree from clustering", VisualConstants.ICON_TYPE_TEXT, "tre");
		manager.plugContentType("application/cel", true, false, "Affymetrix CEL", VisualConstants.ICON_TYPE_RAWDATA, "cel");
	}

	public void plugFeatures(DataManager manager) {
		manager.plugFeatureFactory("/normalised-expression", new NormalisedExpressionProvider());
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
		return SERVER_MODULE_NAME;
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

	@Override
	public boolean isImportToolSupported() {
		return true;
	}

	@Override
	public boolean isWorkflowCompatible(DataBean data) {
		return ChipsterInputTypes.GENE_EXPRS.isTypeOf(data);
	}

	@Override
	public VisualisationMethod[] getVisualisationMethods() {
		return new VisualisationMethod[] {
				VisualisationMethods.ARRAY_LAYOUT,
				VisualisationMethods.HISTOGRAM,
				VisualisationMethods.SCATTERPLOT,
				VisualisationMethods.SCATTERPLOT3D,
				VisualisationMethods.SCATTERPLOT3DPCA,
				VisualisationMethods.VOLCANOPLOT,
				VisualisationMethods.SOM,
				VisualisationMethods.HIERARCHICAL,
				VisualisationMethods.EXPRESSION_PROFILE,
				VisualisationMethods.CLUSTERED_PROFILES,
				VisualisationMethods.VENN_DIAGRAM,
				VisualisationMethods.GBROWSER
		};
	}

	@Override
	public URL getExampleSessionUrl() throws MalformedURLException {
		return new URL(EXAMPLE_SESSION_URL);
	}

	@Override
	public String[][] getRepositoryWorkflows() {
		return new String[][] { 
				{ "Gene expression analysis", "/gene-expression-analysis.bsh" }, 
				{ "Protein expression analysis", "/protein-expression-analysis.bsh" },
				{ "miRNA expression analysis", "/mirna-expression-analysis.bsh" }
		};
	}

}
