package fi.csc.microarray.module.chipster;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.jms.JMSException;
import javax.swing.AbstractAction;
import javax.swing.Icon;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.QuickLinkPanel;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.dialog.TaskImportDialog;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationUtilities;
import fi.csc.microarray.client.visualisation.methods.ArrayLayout;
import fi.csc.microarray.client.visualisation.methods.ChipsterGBrowserVisualisation;
import fi.csc.microarray.client.visualisation.methods.ClusteredProfiles;
import fi.csc.microarray.client.visualisation.methods.ExpressionProfile;
import fi.csc.microarray.client.visualisation.methods.Heatmap;
import fi.csc.microarray.client.visualisation.methods.HierarchicalClustering;
import fi.csc.microarray.client.visualisation.methods.Histogram;
import fi.csc.microarray.client.visualisation.methods.PhenodataEditor;
import fi.csc.microarray.client.visualisation.methods.SOM;
import fi.csc.microarray.client.visualisation.methods.SamBamViewer;
import fi.csc.microarray.client.visualisation.methods.Scatterplot;
import fi.csc.microarray.client.visualisation.methods.Scatterplot3DPCA;
import fi.csc.microarray.client.visualisation.methods.VennDiagram;
import fi.csc.microarray.client.visualisation.methods.Volcanoplot;
import fi.csc.microarray.client.visualisation.methods.threed.Scatterplot3D;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.TypeTag;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.databeans.features.bio.IdentifierProvider;
import fi.csc.microarray.databeans.features.bio.NormalisedExpressionProvider;
import fi.csc.microarray.databeans.features.stat.HierarchicalClusterProvider;
import fi.csc.microarray.databeans.features.stat.SomClusterProvider;
import fi.csc.microarray.databeans.features.table.EditableTable;
import fi.csc.microarray.databeans.features.table.TableBeanEditor;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.LinkUtil;
import fi.csc.microarray.util.Strings;

public class MicroarrayModule implements Module {
	
	public static class TypeTags {
		public static final TypeTag PHENODATA  = new TypeTag("phenodata", "Chipster compatible phenodata");
		public static final TypeTag RAW_AFFYMETRIX_EXPRESSION_VALUES  = new TypeTag("raw-arrymetrix-expression-values", "must be in CEL format (text or binary)");
		public static final TypeTag RAW_EXPRESSION_VALUES  = new TypeTag("raw-expression-values", "");
		public static final TypeTag NORMALISED_EXPRESSION_VALUES = new TypeTag("normalised-expression-values", "must have columns following name pattern \"chip.*\"");
		public static final TypeTag GENENAMES = new TypeTag("genenames", "must have column \" \" or \"identifier\"");
		public static final TypeTag SIGNIFICANT_EXPRESSION_FOLD_CHANGES = new TypeTag("significant-expression-fold-changes", "must have columns following name patterns \"FC.*\" and \"p.*\"");
		public static final TypeTag EXPRESSION_PRIMARY_COMPONENTS_GENEWISE = new TypeTag("expression-primary-components-genewise", "");
		public static final TypeTag EXPRESSION_PRIMARY_COMPONENTS_CHIPWISE = new TypeTag("expression-primary-components-chipwise", "");
		public static final TypeTag ORDERED_GENOMIC_ENTITIES = new TypeTag("ordered-genomic-entities", "");
		public static final TypeTag CLUSTERED_EXPRESSION_VALUES = new TypeTag("clustered-expression-values", "must have column \"cluster\"");
		public static final TypeTag SOM_CLUSTERED_EXPRESSION_VALUES = new TypeTag("som-clustered-expression-values", "must have columns \"colours\", \"distance2first\", \"cluster\", \"griddim\"");
		public static final TypeTag BAM_FILE  = new TypeTag("bam-file", "");
		public static final TypeTag FASTA_FILE  = new TypeTag("fasta-file", "");
		public static final TypeTag GTF_FILE  = new TypeTag("gtf-file", "");
		public static final TypeTag TABLE_WITH_HASH_HEADER = new TypeTag("table-with-hash-header", "header rows start with #");
		public static final TypeTag TABLE_WITH_DOUBLE_HASH_HEADER = new TypeTag("table-with-double-hash-header", "header rows start with ##");
		public static final TypeTag CHROMOSOME_IN_FIRST_TABLE_COLUMN = new TypeTag("chromosome-in-first-table-column", "first column of table is chromosome");
		public static final TypeTag CHROMOSOME_IN_SECOND_TABLE_COLUMN = new TypeTag("chromosome-in-second-table-column", "second column of table is chromosome");
		public static final TypeTag START_POSITION_IN_SECOND_TABLE_COLUMN = new TypeTag("start-position-in-second-table-column", "second column of table is start position");
		public static final TypeTag START_POSITION_IN_THIRD_TABLE_COLUMN = new TypeTag("start-position-in-third-table-column", "third column of table is start position");
		public static final TypeTag START_POSITION_IN_FOURTH_TABLE_COLUMN = new TypeTag("start-position-in-fourth-table-column", "fourth column of table is start position");
		public static final TypeTag END_POSITION_IN_THIRD_TABLE_COLUMN = new TypeTag("end-position-in-third-table-column", "third column of table is end position");
		public static final TypeTag END_POSITION_IN_FOURTH_TABLE_COLUMN = new TypeTag("end-position-in-fourth-table-column", "fourth column of table is end position");
		public static final TypeTag END_POSITION_IN_FIFTH_TABLE_COLUMN = new TypeTag("end-position-in-fifth-table-column", "fifth column of table is end position");
		public static final TypeTag CNA = new TypeTag("cna-data", "data file produced by cna tools");
		public static final TypeTag MOTHUR_OLIGOS = new TypeTag("Mothur oligos data", "Mothur oligos data");
		public static final TypeTag MOTHUR_NAMES = new TypeTag("Mothur names data", "Mothur names data");
		public static final TypeTag MOTHUR_GROUPS = new TypeTag("Mothur groups data", "Mothur groups data");
	}
	
	public static class VisualisationMethods {
		public static VisualisationMethod ARRAY_LAYOUT = new VisualisationMethod("Array layout", ArrayLayout.class, VisualConstants.ARRAY_MENUICON, 100, 0.0009);
		public static VisualisationMethod HISTOGRAM = new VisualisationMethod("Histogram", Histogram.class, VisualConstants.HISTOGRAM_MENUICON, 70, 0.024);
		public static VisualisationMethod HEATMAP = new VisualisationMethod("Heatmap", Heatmap.class, VisualConstants.HEATMAP_MENUICON, 100, 0.0009);
		public static VisualisationMethod SCATTERPLOT = new VisualisationMethod("Scatterplot", Scatterplot.class, VisualConstants.SCATTER_MENUICON, 79, 0.039);
		public static VisualisationMethod SCATTERPLOT3D = new VisualisationMethod("3D Scatterplot", Scatterplot3D.class, VisualConstants.SCATTER3D_MENUICON, 78, 0.082);
		public static VisualisationMethod SCATTERPLOT3DPCA = new VisualisationMethod("3D Scatterplot for PCA", Scatterplot3DPCA.class, VisualConstants.SCATTER3DPCA_MENUICON, 77, 0.082);
		public static VisualisationMethod VOLCANOPLOT = new VisualisationMethod("Volcano plot", Volcanoplot.class, VisualConstants.VOLCANO_MENUICON, 80, 0.039); 
		public static VisualisationMethod SOM = new VisualisationMethod("SOM", SOM.class, VisualConstants.SOM_MENUICON, 59, 0.034);
		public static VisualisationMethod HIERARCHICAL = new VisualisationMethod("Hierarchical clustering", HierarchicalClustering.class, VisualConstants.HC_MENUICON, 60, 0.09);
		public static VisualisationMethod EXPRESSION_PROFILE = new VisualisationMethod("Expression profile", ExpressionProfile.class, VisualConstants.PROFILE_MENUICON, 90, 0.1);
		public static VisualisationMethod CLUSTERED_PROFILES = new VisualisationMethod("Clustered profiles", ClusteredProfiles.class, VisualConstants.PROFILES_MENUICON, 89, 0.087);
		public static VisualisationMethod VENN_DIAGRAM = new VisualisationMethod("Venn-diagram", VennDiagram.class, VisualConstants.VENN_MENUICON, 50, -1);
		public static VisualisationMethod GBROWSER = new VisualisationMethod("Genome browser", ChipsterGBrowserVisualisation.class, VisualConstants.GB_MENUICON, 90, -1, "genomeBrowser.html");
		public static VisualisationMethod SAMBAM_VIEWER = new VisualisationMethod("BAM viewer", SamBamViewer.class, VisualConstants.TEXT_MENUICON, 100, -1);
		public static VisualisationMethod PHENODATA = new VisualisationMethod("Phenodata editor", PhenodataEditor.class, VisualConstants.PHENODATA_MENUICON, 105, 0, "visualisation-phenodata.html");
	}
	
	public static final String SERVER_MODULE_NAME = "microarray";

	public static final String ANNOTATION_ID = "annotate-genelist2html.R";

	public static final String IMPORT_FROM_ARRAYEXPRESS_ID = "import-ArrayExpress.R";
	public static final String IMPORT_FROM_GEO_ID = "import-soft2.R";

	public void plugContentTypes(DataManager manager) {
		manager.plugContentType("application/x-treeview", true, false, "Newick formatted tree from clustering", VisualConstants.ICON_TYPE_TEXT, "tre");
		manager.plugContentType("application/cel", true, false, "Affymetrix CEL", VisualConstants.ICON_TYPE_RAWDATA, "cel");
		manager.plugContentType("text/bed", true, false, "BED file", VisualConstants.ICON_TYPE_TEXT, "bed");
		manager.plugContentType("text/gtf", true, false, "Gene Transfer Format file", VisualConstants.ICON_TYPE_TEXT, "gtf", "gff", "gff2", "gff3");
		manager.plugContentType("chemical/x-fasta", true, false, "FASTA", VisualConstants.ICON_TYPE_TEXT, "fasta", "fa", "fna", "fsa", "mpfa");
		manager.plugContentType("text/fastq", true, false, "FASTQ", VisualConstants.ICON_TYPE_TEXT, "fastq", "fq");
		manager.plugContentType("application/gzip", true, true, "Gzip file", VisualConstants.ICON_TYPE_BINARY, "gz");
		manager.plugContentType("text/vcf", true, false, "Variant Call Format", VisualConstants.ICON_TYPE_TEXT, "vcf");
		manager.plugContentType("application/bam", true, false, "Binary sequence Alignment/Map format", VisualConstants.ICON_TYPE_TEXT, "bam");
		manager.plugContentType("text/qual", true, false, "Quality file", VisualConstants.ICON_TYPE_TEXT, "qual");
		manager.plugContentType("text/mothur-oligos", true, false, "Mothur oligos file", VisualConstants.ICON_TYPE_TEXT, "oligos");
		manager.plugContentType("text/mothur-names", true, false, "Mothur names file", VisualConstants.ICON_TYPE_TEXT, "names");
		manager.plugContentType("text/mothur-groups", true, false, "Mothur groups file", VisualConstants.ICON_TYPE_TEXT, "groups");
		manager.plugContentType("text/sff", true, false, "sff file", VisualConstants.ICON_TYPE_TEXT, "sff");
	}
	

	public void plugFeatures(DataManager manager) {
		manager.plugFeatureFactory("/normalised-expression", new NormalisedExpressionProvider());
		manager.plugFeatureFactory("/identifier", new IdentifierProvider());;
		manager.plugFeatureFactory("/clusters/som", new SomClusterProvider());
		manager.plugFeatureFactory("/clusters/hierarchical", new HierarchicalClusterProvider());
	}

	public void plugModifiers(DataManager manager) {
		// nothing to plug
	}

	@Override
	public String[] getServerModuleNames() {
		return new String[] { "microarray", "ngs", "misc" };
	}

	@Override
	public String getModuleLongName(String moduleName) {
		if ("microarray".equals(moduleName)) {
			return "Microarrays";
		} else if ("ngs".equals(moduleName)) {
			return "NGS";
		} else if ("misc".equals(moduleName)) {
			return "Misc";
		} else {
			return moduleName;
		}
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

	@SuppressWarnings("serial")
	@Override
	public void addImportLinks(QuickLinkPanel quickLinkPanel, List<JXHyperlink> importLinks) {

		importLinks.add(LinkUtil.createLink("Import from ArrayExpress", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent e) {
				doImportFromArrayExpress();
			}
		}));

		importLinks.add(LinkUtil.createLink("Import from GEO", new AbstractAction() {
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
			new TaskImportDialog(application, "Import data from the GEO", null, importOperation);
			
		} catch (Exception me) {
			Session.getSession().getApplication().reportException(me);
		}
	}

	private void doImportFromArrayExpress() {
		try {
			ClientApplication application = Session.getSession().getApplication();
			Operation importOperation = new Operation(application.getOperationDefinition(MicroarrayModule.IMPORT_FROM_ARRAYEXPRESS_ID), new DataBean[] {});
			new TaskImportDialog(application, "Import data from the ArrayExpress", null, importOperation);
			
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
		// TODO replace inclusive check with exclusive: check for raw expression data and other illegal types
//		return ChipsterInputTypes.GENE_EXPRS.isTypeOf(data);
		return true;
	}

	@Override
	public VisualisationMethod[] getVisualisationMethods() {
		return new VisualisationMethod[] {
				VisualisationMethods.PHENODATA,
				VisualisationMethods.ARRAY_LAYOUT,
				VisualisationMethods.HEATMAP,
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
				VisualisationMethods.SAMBAM_VIEWER,
				VisualisationMethods.GBROWSER
		};
	}

	@Override
	public List<DbSession> getExampleSessions(boolean isStandalone) throws JMSException, Exception {
		
		 List<DbSession> sessions = Session.getSession().getServiceAccessor().getFileBrokerClient().listPublicRemoteSessions();
		 
		return sessions;
	}

	@Override
	public String[][] getRepositoryWorkflows() {
		return new String[][] { 
				{ "Gene expression analysis", "/gene-expression-analysis.bsh" }, 
				{ "Protein expression analysis", "/protein-expression-analysis.bsh" },
				{ "miRNA expression analysis", "/mirna-expression-analysis.bsh" }
		};
	}
	
	
	/**
	 * Selects the proper source dataset ie. the dataset that is not a hidden phenodata dataset.
	 */
	public static DataBean getProperSource(DataBean dataBean) {
		
		if (dataBean == null || dataBean.getLinkTargets(Link.DERIVATION).size() == 0) {
			return null;
			
		} else if (dataBean.getLinkTargets(Link.DERIVATION).size() == 1) {
			return dataBean.getLinkTargets(Link.DERIVATION).iterator().next();
			
		} else {
			LinkedList<DataBean> sourceCollector = new LinkedList<DataBean>();
			for (DataBean source : dataBean.getLinkTargets(Link.DERIVATION)) {
				if (source.queryFeatures("/phenodata").exists()) {
					sourceCollector.add(source);
				}
			}
			if (sourceCollector.size() == 0 || sourceCollector.size() > 1) {
				return null; // no definite source was found
			}
			
			return sourceCollector.getFirst(); // return the one and only proper source
		}
	}
	
	/**
	 * Gets the operation history of this dataset as a chronological list
	 * (implemented as vector) of DataSetHistoryWrapper objects, which
	 * practically are DataSets with a slightly altered toString output.
	 * The list will contain all the datasets on the workflow, starting
	 * from raw data and ending at current dataset.
	 * 
	 * @return The chronologically ascending list of dataset history wrappers.
	 */
	public static DataBean[] getSourcePath(DataBean dataBean) {
		LinkedList<DataBean> list = new LinkedList<DataBean>();
		if (getProperSource(dataBean) != null) {
			list.addAll(Arrays.asList(getSourcePath(getProperSource(dataBean))));
		}
		list.add(dataBean);
		return list.toArray(new DataBean[0]);
	}	


	@Override
	public boolean isMetadata(DataBean data) {
		// FIXME how this should actually work in the new type system?
		return data.queryFeatures("/phenodata").exists();
	}
	
	@Override
	public void preProcessInputMetadata(Operation oper, DataBean metadataInput) throws MicroarrayException, IOException {
		
		// R can't handle spaces in column titles. Terminate job if there are any white space characters.
		
		Table table = metadataInput.queryFeatures("/column/*").asTable();
		
		Pattern pattern = Pattern.compile("\\s"); // whitespace characters: space, tab and new line
		for (String column : table.getColumnNames()) {						
			Matcher matcher = pattern.matcher(column);
			if (matcher.find()) { 
				throw new MicroarrayException("Phenodata columnn titles must not contain white space characters (space, tab and new line)");
			}
		}
	}

	@Override
	public void postProcessOutputMetadata(Operation oper, DataBean metadataOutput) throws MicroarrayException, IOException {
		
		// FIXME how this should actually work in the new type system?
		
		// if original names are not already contained in the phenodata
		if (!metadataOutput.queryFeatures("/column/" + PhenodataEditor.PHENODATA_NAME_COLUMN).exists()) {
			// augment phenodata with original dataset names (using parameter bindings)
			HashSet<String> insertedNames = new HashSet<String>();
			TableBeanEditor tableEditor = new TableBeanEditor(metadataOutput);
			EditableTable editableTable = tableEditor.getEditable();
			LinkedList<String> newColumn = new LinkedList<String>();
			newColumn.addAll(Arrays.asList(Strings.repeatToArray("", editableTable.getRowCount())));
			editableTable.addColumn(PhenodataEditor.PHENODATA_NAME_COLUMN, 1, newColumn); // add after sample column 
			for (int ri = 0; ri < editableTable.getRowCount(); ri++) {
				String sample = editableTable.getValue(PhenodataEditor.PHENODATA_SAMPLE_COLUMN, ri);
				boolean correctRowFound = false;
				String originalName = null;
				for (DataBinding binding : oper.getBindings()) {
					if (binding.getName().equals(sample)) {
						DataBean ancestor = binding.getData().getUniqueAncestorRecursively(binding.getData());
						if (!ancestor.equals(binding.getData())) {
							originalName = binding.getData().getName() + " ( " + ancestor.getName() + " )";
						} else {
							originalName = binding.getData().getName();
						}
						correctRowFound = true;
						break;
					}
				}
				if (!correctRowFound) {
					originalName = sample; // just duplicate the sample name if proper is not found
				}

				// check that original names are unique
				if (insertedNames.contains(originalName)) {
					final String separator = "/";
					int i = 2;
					while (insertedNames.contains(originalName + separator + i)) {
						i++;
					}
					originalName = originalName + separator + i;
				}

				editableTable.setValue(PhenodataEditor.PHENODATA_NAME_COLUMN, ri, originalName);
				insertedNames.add(originalName);

			}
			tableEditor.write();
		}

		// if chip descriptions (visualisation view names) are not already contained in the phenodata
		if (!metadataOutput.queryFeatures("/column/" + PhenodataEditor.PHENODATA_DESCRIPTION_COLUMN).exists()) {
			// copy original dataset names
			TableBeanEditor tableEditor = new TableBeanEditor(metadataOutput);
			EditableTable editableTable = tableEditor.getEditable();
			LinkedList<String> newColumn = new LinkedList<String>();
			newColumn.addAll(Arrays.asList(Strings.repeatToArray("", editableTable.getRowCount())));
			editableTable.addColumn(PhenodataEditor.PHENODATA_DESCRIPTION_COLUMN, newColumn); 
			for (int ri = 0; ri < editableTable.getRowCount(); ri++) {
				String sample = editableTable.getValue(PhenodataEditor.PHENODATA_NAME_COLUMN, ri);										
				editableTable.setValue(PhenodataEditor.PHENODATA_DESCRIPTION_COLUMN, ri, sample);
			}
			tableEditor.write();
		}
	}

	@Override
	public String getShortDataName(String categoryName) {
		return BasicModule.shortenDataName(categoryName);
	}

	@Override
	public boolean countOperationResults() {
		return true;
	}

	/**
	 * Generates nice context link panel for quickly using genome browser. If not in standalone
	 * mode, null is returned. 
	 */
	@SuppressWarnings("serial")
	@Override
	public JPanel getContextLinkPanel(int selectedDataCount) {

		final ClientApplication application = Session.getSession().getApplication();

		if (!application.isStandalone()) {
			return null;
		}

		JPanel contentPanel = new JPanel();

		// Initialise context link panel
		contentPanel.setBackground(Color.WHITE);
		contentPanel.setLayout(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();
		c.gridx = 0;
		c.gridy = 0;
		c.anchor = GridBagConstraints.NORTHWEST;

		int topMargin = 15;
		int leftMargin = 30;
		c.insets.set(topMargin, leftMargin, 0, 0);
		contentPanel.add(new JLabel(VisualConstants.QUICKLINK_ICON), c);
		JXHyperlink link;
		boolean currentSelectionVisualisable = false;
		try {
			Visualisation visualisation = MicroarrayModule.VisualisationMethods.GBROWSER.getVisualiser(null);
			if (visualisation != null) {
				List<DataBean> selection = application.getSelectionManager().getSelectedDataBeans();
				currentSelectionVisualisable = visualisation.canVisualise(selection);
			}
		} catch (Exception e2) {
			// ignore
		}

		if (selectedDataCount > 0 && currentSelectionVisualisable) {
			// Selection can be visualised directly
			link = new JXHyperlink(new AbstractAction() {
				@Override
				public void actionPerformed(ActionEvent e) {
					application.setVisualisationMethod(MicroarrayModule.VisualisationMethods.GBROWSER, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
				}
			});
			link.setText("Open genome browser");


		} else if (selectedDataCount < application.getDataManager().databeans().size()) {
			// Not all datasets are selected, they might be visualisable when selected, so suggest it
			link = new JXHyperlink(new AbstractAction() {
				@Override
				public void actionPerformed(ActionEvent e) {

					// Select all datasets (creates a bunch of events)
					application.selectAllItems();

					// Check if selection can be visualised
					boolean canVisualise = false;
					try {
						canVisualise = MicroarrayModule.VisualisationMethods.GBROWSER.getVisualiser(null).canVisualise(application.getSelectionManager().getSelectedDataBeans());
					} catch (Exception e1) {
						// ignore
					}

					if (canVisualise) {
						// Use invokeLater so that visualisation system can process events from selectAllItems 
						// before this. Otherwise the system is not synchronized and the visualisation panel
						// does not update.
						SwingUtilities.invokeLater(new Runnable() {
							@Override
							public void run() {
								application.setVisualisationMethod(MicroarrayModule.VisualisationMethods.GBROWSER, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
							}
						});
					} else {
						application.showDialog("Cannot open genome browser", "Currect dataset selection cannot be visualised in genome browser", "Unselect improper datasets to open genome browser.", Severity.INFO, true);
					}
				}
			});
			link.setText("Select all and open genome browser");

		} else {
			// No links available, bail out
			return null;
		}

		c.insets.set(topMargin, 5, 0, 0);
		c.gridx++;
		contentPanel.add(link, c);

		return contentPanel;
	}

	@Override
	public boolean notesVisibleAtStartup() {
		return true;
	}

	@Override
	public String getDisplayName() {
		return "Chipster";
	}

	@Override
	public String getManualHome() {
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		return configuration.getString("client", "manual-root");
	}

	@Override
	public void addSpeadsheetMenuItems(JPopupMenu spreadsheetPopupMenu, final VisualisationFrame visualisationFrame) {

		JMenuItem annotateMenuItem = new JMenuItem();
		annotateMenuItem.setText("Create dataset and annotate with Bioconductor");
		annotateMenuItem.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent e) {
				VisualisationUtilities.annotateBySelection(visualisationFrame.getDatas(), ANNOTATION_ID);
			}
		});
		spreadsheetPopupMenu.add(annotateMenuItem);
	}

	@Override
	public List<Boolean> flagLinkableColumns(Table columns, DataBean data) {
		LinkedList<Boolean> flags = new LinkedList<Boolean>();
		
		flags.add(data.hasTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN)); //0
		flags.add(data.hasTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_SECOND_TABLE_COLUMN) ||
				data.hasTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_SECOND_TABLE_COLUMN));  //1
		flags.add(data.hasTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_THIRD_TABLE_COLUMN) ||
				data.hasTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_THIRD_TABLE_COLUMN)); //2
		flags.add(data.hasTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FOURTH_TABLE_COLUMN)||
				data.hasTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_FOURTH_TABLE_COLUMN)); //3
		flags.add(data.hasTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FIFTH_TABLE_COLUMN)); //4
		
		for (int i = 4; i < columns.getColumnCount(); i++) {
			flags.add(false);
		}
		return flags;
	}

	@Override
	public IntegratedEntity createLinkableEntity(Table columns, DataBean data) {
		
		IntegratedEntity entity = new IntegratedEntity();
		
		final String CHR_KEY = "chromosome";
		final String START_KEY = "start";
		final String END_KEY = "end";
		
		
		//Chromosome
		int chrColumn = -1;
		
		if (data.hasTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN)) {
			chrColumn = 0;
		}
		if (data.hasTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_SECOND_TABLE_COLUMN)) {
			chrColumn = 1;
		}
		
		if (chrColumn != -1) {
			entity.put(CHR_KEY, columns.getStringValue(columns.getColumnNames()[chrColumn]));
		}
		
		//Start
		int startColumn = -1;
		
		if (data.hasTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_SECOND_TABLE_COLUMN)) {
			startColumn = 1;
		}
		
		if (data.hasTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_THIRD_TABLE_COLUMN)) {
			startColumn = 2;
		}

		if (data.hasTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_FOURTH_TABLE_COLUMN)) {
			startColumn = 3;
		}

		if (startColumn != -1) {
			entity.put(START_KEY, columns.getStringValue(columns.getColumnNames()[startColumn]));
		}
		
		//End
		int endColumn = -1;
		
		if (data.hasTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_THIRD_TABLE_COLUMN)) {
			endColumn = 2;
		}

		if (data.hasTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FOURTH_TABLE_COLUMN)) {
			endColumn = 3;
		}
		
		if (data.hasTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FIFTH_TABLE_COLUMN)) {
			endColumn = 4;
		}
		
		if (endColumn != -1) {
			entity.put(END_KEY, columns.getStringValue(columns.getColumnNames()[endColumn]));
		}		

		return entity;
	}

	@Override
	public void addTypeTags(DataBean data) throws MicroarrayException, IOException {

		if (data.isContentTypeCompatitible("application/cel")) {
			data.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
		}
		
		if (data.isContentTypeCompatitible("text/tab", "application/cel", "text/csv")) {
			data.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
		}
		
		if (data.isContentTypeCompatitible("text/tab")) {
			BufferedReader in = null;
			try {
				in = new BufferedReader(new InputStreamReader(Session.getSession().getDataManager().getContentStream(data, DataNotAvailableHandling.EXCEPTION_ON_NA)));
				String headerLine = in.readLine();
				String contentLine = in.readLine();
				
				if (headerLine != null && contentLine != null) {						
					
					List<String> chrColumns = Arrays.asList(new String[] { "chr", "SEQNAMES", "chromosome"});
					List<String> startColumns = Arrays.asList(new String[] { "start", "START"});
					List<String> endColumns = Arrays.asList(new String[] { "end", "END"});											

					String[] split = headerLine.split("\t");
					
					int shiftHeader = 0;
					
					//TableColumnProvider has this check already, but can we use Feature api here?
					if (split.length == contentLine.split("\t").length - 1) {
						shiftHeader = 1;
					}
					
					for (int i = 0; i < split.length; i++) {
						String column = split[i];
					
						if (chrColumns.contains(column)) {
							switch(i + shiftHeader) {
							case 0:
								data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN);
								break;
							case 1:
								data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_SECOND_TABLE_COLUMN);
								break;
							}
						}
						if (startColumns.contains(column)) {
							switch(i + shiftHeader) {
							case 1:
								data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_SECOND_TABLE_COLUMN);
								break;
							case 2:									
								data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_THIRD_TABLE_COLUMN);
								break;
							}
						}
						if (endColumns.contains(column)) {
							switch(i + shiftHeader) {
							case 2:
								data.addTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_THIRD_TABLE_COLUMN);
								break;
							case 3:
								data.addTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FOURTH_TABLE_COLUMN);
								break;
							}
						}
						
						if (column.contains("loss.freq") ||
							column.contains("gain.freq") ||
							column.startsWith("flag.") || 
							column.startsWith("segmented.")) {
							
							if (!data.hasTypeTag(MicroarrayModule.TypeTags.CNA)) {
								data.addTypeTag(MicroarrayModule.TypeTags.CNA);
							}
						}							
					}						 
				}
				
			} catch (IOException e) {
				throw new RuntimeException(e);
				
			} finally {
				IOUtils.closeIfPossible(in);
			}				
		}

		if (data.isContentTypeCompatitible("text/gtf")) {
			data.addTypeTag(MicroarrayModule.TypeTags.GTF_FILE);
			data.addTypeTag(BasicModule.TypeTags.TABLE_WITHOUT_COLUMN_NAMES);
			data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN);
			data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_FOURTH_TABLE_COLUMN);				
			data.addTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FIFTH_TABLE_COLUMN);				
			addTypeTagIfHashHeader(data);
		}
		
		if (data.isContentTypeCompatitible("text/vcf")) {
			data.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
			addTypeTagIfDoubleHashHeader(data);
			data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN);
			data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_SECOND_TABLE_COLUMN);
		}

		if (data.isContentTypeCompatitible("text/bed")) {
			data.addTypeTag(BasicModule.TypeTags.TABLE_WITHOUT_COLUMN_NAMES);
			
			data.addTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN);
			data.addTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_SECOND_TABLE_COLUMN);				
			data.addTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_THIRD_TABLE_COLUMN);

			if (readFirstLine(data).startsWith("track")) {
				data.addTypeTag(BasicModule.TypeTags.TABLE_WITH_HEADER_ROW);
			}
		}

		if (data.isContentTypeCompatitible("application/bam")) {
			data.addTypeTag(MicroarrayModule.TypeTags.BAM_FILE);
		}
		
		if (data.isContentTypeCompatitible("chemical/x-fasta")) {
            data.addTypeTag(MicroarrayModule.TypeTags.FASTA_FILE);
		}
		
        // mothur
        if (data.isContentTypeCompatitible("text/mothur-oligos")) {
                data.addTypeTag(MicroarrayModule.TypeTags.MOTHUR_OLIGOS);
        }
        if (data.isContentTypeCompatitible("text/mothur-names")) {
                data.addTypeTag(MicroarrayModule.TypeTags.MOTHUR_NAMES);
        }
        if (data.isContentTypeCompatitible("text/mothur-groups")) {
                data.addTypeTag(MicroarrayModule.TypeTags.MOTHUR_GROUPS);
        }

		// Rest of the tags are set only when this module is primary
		
		if (!(Session.getSession().getPrimaryModule() instanceof MicroarrayModule)) {
			return;
		}

		Table chips = data.queryFeatures("/column/chip.*").asTable();

		if (data.isContentTypeCompatitible("application/cel")) {
			data.addTypeTag(MicroarrayModule.TypeTags.RAW_AFFYMETRIX_EXPRESSION_VALUES);

		} else if (data.queryFeatures("/column/sample").exists() && !data.queryFeatures("/phenodata").exists()) {
			data.addTypeTag(MicroarrayModule.TypeTags.RAW_EXPRESSION_VALUES);

		} else if (chips != null && chips.getColumnCount() > 0) {
			data.addTypeTag(MicroarrayModule.TypeTags.NORMALISED_EXPRESSION_VALUES);
		} 

		if (data.queryFeatures("/identifier").exists()) {
			data.addTypeTag(MicroarrayModule.TypeTags.GENENAMES);
		} 


		// Tag additional typing information
		if (data.queryFeatures("/phenodata").exists()) {
			data.addTypeTag(MicroarrayModule.TypeTags.PHENODATA);
		}

		if (data.queryFeatures("/column/p.*").exists() && data.queryFeatures("/column/FC*").exists()) {
			data.addTypeTag(MicroarrayModule.TypeTags.SIGNIFICANT_EXPRESSION_FOLD_CHANGES);
		}

		boolean isChipwise = false;
		ParameterRecord pcaOn = data.getOperationRecord().getParameter("do.pca.on");
		if (pcaOn != null) {
			String pcaOnValue = pcaOn.getValue();
			if (pcaOnValue != null && pcaOnValue.equals("chips")) {
				isChipwise = true;
			}
		}
		if (data.getOperationRecord().getNameID().getID().equals("ordination-pca.R") && isChipwise) {
			data.addTypeTag(MicroarrayModule.TypeTags.EXPRESSION_PRIMARY_COMPONENTS_CHIPWISE);
		}

		if (chips != null && chips.getColumnNames().length > 1 && data.queryFeatures("/column/cluster").exists()) {
			data.addTypeTag(MicroarrayModule.TypeTags.CLUSTERED_EXPRESSION_VALUES);
		}

		if (data.queryFeatures("/clusters/som").exists()) {
			data.addTypeTag(MicroarrayModule.TypeTags.SOM_CLUSTERED_EXPRESSION_VALUES);
		}

		// Finally, set NGS related tags
		if (data.isContentTypeCompatitible("text/bed") 
				|| (data.isContentTypeCompatitible("application/octet-stream")) && (data.getName().contains(".bam-summary")) 
				|| (data.isContentTypeCompatitible("application/octet-stream")) && (data.getName().endsWith(".bam") || data.getName().endsWith(".sam"))
				|| (data.isContentTypeCompatitible("application/octet-stream")) && (data.getName().endsWith(".bai"))) {

			data.addTypeTag(MicroarrayModule.TypeTags.ORDERED_GENOMIC_ENTITIES);
			
		} else if (data.isContentTypeCompatitible("text/tab")) {
			
			// require .tsv to have columns for genomic coordinates
			String line = readFirstLine(data); 
			if (line != null && line.contains("chr") && line.contains("start") && line.contains("end")) { 
				data.addTypeTag(MicroarrayModule.TypeTags.ORDERED_GENOMIC_ENTITIES);
			}
			
		}

	}
	
	private String readFirstLine(DataBean data) {

		BufferedReader reader = null;
		try {
			InputStream stream = Session.getSession().getDataManager().getContentStream(data, DataNotAvailableHandling.NULL_ON_NA);
			if (stream != null) {
				reader = new BufferedReader(new InputStreamReader(stream));
				return reader.readLine();
			} else {
				return "";
			}
			
		} catch (IOException e) {
			throw new RuntimeException(e);
		} finally {
			IOUtils.closeIfPossible(reader);
		}
	}

	@Override
	public Icon getIconFor(DataBean data) {
		if (data.hasTypeTag(MicroarrayModule.TypeTags.PHENODATA)) {
			return VisualConstants.ICON_TYPE_PHENODATA;
		} else {
			return data.getContentType().getIcon();
		}
	}

	private void addTypeTagIfHashHeader(DataBean data) {
		String line = readFirstLine(data);
		if (line != null && line.startsWith("#")) {
			data.addTypeTag(MicroarrayModule.TypeTags.TABLE_WITH_HASH_HEADER);
		}
	}

	private void addTypeTagIfDoubleHashHeader(DataBean data) {
		String line = readFirstLine(data);
		if (line != null && line.startsWith("##")) {
			data.addTypeTag(MicroarrayModule.TypeTags.TABLE_WITH_DOUBLE_HASH_HEADER);
		}
	}

}
