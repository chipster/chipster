package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.CardLayout;
import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;

import org.jdesktop.swingx.JXHyperlink;
import org.jfree.chart.JFreeChart;

import fi.csc.chipster.tools.gbrowser.SamBamUtils;
import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.LinkUtil;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.selection.PointSelectionEvent;
import fi.csc.microarray.client.visualisation.NonScalableChartPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomePlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.BedTabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GeneSearchHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GtfTabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixSummaryHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParserWithCoordinateConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.HeaderTsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.VcfParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.AnnotationType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.gbrowser.index.GeneIndexActions;
import fi.csc.microarray.util.BrowserLauncher;
import fi.csc.microarray.util.IOUtils;

/**
 * Facade class that hides genome browser internals and exposes an API that is compatible 
 * with Chipster visualization system. 
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * @see GenomePlot
 */
public class GenomeBrowser extends Visualisation implements ActionListener,
RegionListener, ComponentListener, PropertyChangeListener {


	private static final long DEFAULT_VIEWSIZE = 100000;
	private static final long DEFAULT_LOCATION = 1000000;
	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";
	private static final String COVERAGE_NONE = "none";
	private static final String COVERAGE_TOTAL = "total";
	private static final String COVERAGE_STRAND = "strand-specific";

	private static class Interpretation {

		public TrackType type;
		public List<DataBean> summaryDatas = new LinkedList<DataBean>();
		public DataBean primaryData;
		public DataBean indexData;

		public Interpretation(TrackType type, DataBean primaryData) {
			this.type = type;
			this.primaryData = primaryData;
		}

	}

	private static enum TrackType {
		CYTOBANDS(false), 
		GENES(false), 
		TRANSCRIPTS(true), 
		REFERENCE(true),
		REGIONS(true),
		REGIONS_WITH_HEADER(true), 
		READS(true),
		HIDDEN(false), 
		VCF(true);

		private boolean isToggleable;

		private TrackType(boolean toggleable) {
			this.isToggleable = toggleable;
		}
	}

	private static class Track {

		Interpretation interpretation;
		JCheckBox checkBox;
		String name;
		TrackGroup trackGroup = null;

		public Track(String name, Interpretation interpretation) {
			this.name = name;
			this.interpretation = interpretation;
		}

		public void setTrackGroup(TrackGroup trackGroup) {
			this.trackGroup = trackGroup;
		}
	}

	private final ClientApplication application = Session.getSession()
			.getApplication();

	private List<Interpretation> interpretations;
	private List<Track> tracks = new LinkedList<Track>();

	private GenomePlot plot;

	private JPanel paramPanel;
	private JPanel settingsPanel = new JPanel();
	private JPanel genomePanel;
	private JPanel locationPanel;
	private JPanel datasetsPanel;
	private JPanel datasetSwitchesPanel;
	private JPanel optionsPanel;
	private JPanel linksPanel;
	private JPanel plotPanel = new JPanel(new CardLayout());

	private JButton goButton = new JButton("Go");

	private JLabel locationLabel = new JLabel("Location (gene or position)");
	private JTextField locationField = new JTextField();

	private JLabel viewsizeLabel = new JLabel("View size");
	private JTextField viewsizeField = new JTextField();

	private JLabel chrLabel = new JLabel("Chromosome");
	private JComboBox chrBox = new JComboBox();

	private JComboBox genomeBox = new JComboBox();

	private AnnotationManager annotationManager;

	private JLabel coverageScaleLabel = new JLabel("Coverage scale");
	private JComboBox coverageScaleBox = new JComboBox();

	private JLabel coverageTypeLabel = new JLabel("Coverage type");
	private JComboBox coverageTypeBox = new JComboBox(); 

	private GeneIndexActions gia;

	private boolean initialised;

	private Map<JCheckBox, String> trackSwitches = new LinkedHashMap<JCheckBox, String>();
	private Set<JCheckBox> datasetSwitches = new HashSet<JCheckBox>();
	private Long lastLocation;
	private Long lastViewsize;
	private JScrollPane verticalScroller;
	private JCheckBox showFullHeightBox;
	private ViewLimiter viewLimiter;
	protected boolean geneSearchDone;
	private JXHyperlink ensemblLink;
	private JXHyperlink ucscLink;


	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);

		// initialize annotations
		this.annotationManager = new AnnotationManager();
		this.annotationManager.initialize();

		trackSwitches.put(new JCheckBox("Reads", true), "Reads");
		trackSwitches.put(new JCheckBox("Highlight SNPs", true), "highlightSNP");
		//		trackSwitches.put(new JCheckBox("Coverage and SNPs", true), "ProfileSNPTrack");
		//		trackSwitches.put(new JCheckBox("Strand-specific coverage", false), "ProfileTrack");

		//		trackSwitches.put(new JCheckBox("Quality coverage", false), "QualityCoverageTrack"); // TODO re-enable quality coverage
		trackSwitches.put(new JCheckBox("Density graph", false), "GelTrack");
		trackSwitches.put(new JCheckBox("Low complexity regions", false), "RepeatMaskerTrack");
		//		trackSwitches.put(new JCheckBox("Known SNP's", false), "changeSNP"); // TODO re-enable dbSNP view		
	}

	@Override
	public JPanel getParameterPanel() {

		if (paramPanel == null) {
			paramPanel = new JPanel();
			paramPanel.setLayout(new GridBagLayout());

			JPanel settings = this.createSettingsPanel();
			JScrollPane settingsScrollPane = new JScrollPane(settings);
			settingsScrollPane.setBorder(BorderFactory.createEmptyBorder());
			settingsScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Settings", settingsScrollPane);
			tabPane.addTab("Legend", new GBrowserLegend());

			GridBagConstraints c = new GridBagConstraints();

			c.gridy = 0;
			c.gridx = 0;
			c.anchor = GridBagConstraints.NORTHWEST;
			c.fill = GridBagConstraints.BOTH;
			c.weighty = 1.0;
			c.weightx = 1.0;
			c.insets.set(5, 0, 0, 0);

			paramPanel.add(tabPane, c);
		}
		
		return paramPanel;
	}

	private void createAvailableTracks() {

		// for now just always add genes and cytobands
		tracks.add(new Track(AnnotationManager.AnnotationType.GTF_TABIX.getId(), new Interpretation(TrackType.GENES, null)));
		tracks.add(new Track(AnnotationManager.AnnotationType.CYTOBANDS.getId(), new Interpretation(TrackType.CYTOBANDS, null)));


		for (int i = 0; i < interpretations.size(); i++) {
			Interpretation interpretation = interpretations.get(i);
			DataBean data = interpretation.primaryData;
			tracks.add(new Track(data.getName(), interpretation));
		}

		// update the dataset switches in the settings panel
		updateDatasetSwitches();
	}

	private JPanel getDatasetsPanel() {
		if (datasetsPanel == null) { 
			datasetsPanel = new JPanel(new GridBagLayout());
			datasetsPanel.setBorder(VisualConstants.createSettingsPanelSubPanelBorder("Datasets"));

			GridBagConstraints dc = new GridBagConstraints();
			dc.gridy = 0;
			dc.gridx = 0;
			dc.anchor = GridBagConstraints.NORTHWEST;
			dc.fill = GridBagConstraints.BOTH;
			dc.weighty = 1.0;
			dc.weightx = 1.0;

			datasetSwitchesPanel = new JPanel();
			datasetSwitchesPanel.setLayout(new BoxLayout(datasetSwitchesPanel, BoxLayout.Y_AXIS));

			datasetsPanel.add(datasetSwitchesPanel, dc);
		}

		return datasetsPanel;
	}

	private void updateDatasetSwitches() {
		datasetSwitches.clear();

		datasetSwitchesPanel.removeAll();

		for (Track track : tracks) {
			JCheckBox box = new JCheckBox(track.name, true);
			box.setToolTipText(track.name);
			box.setEnabled(false);
			track.checkBox = box;
			if (track.interpretation.type.isToggleable) {
				datasetSwitchesPanel.add(box);
				datasetSwitches.add(box);
				box.addActionListener(this);
			}
		}
	}

	public JPanel getOptionsPanel() {
		if (this.optionsPanel == null) {
			optionsPanel = new JPanel(new GridBagLayout());
			optionsPanel.setBorder(VisualConstants.createSettingsPanelSubPanelBorder("Options"));

			GridBagConstraints oc = new GridBagConstraints();
			oc.gridy = 0;
			oc.gridx = 0;
			oc.anchor = GridBagConstraints.PAGE_START;
			oc.fill = GridBagConstraints.BOTH;
			oc.weighty = 1.0;
			oc.weightx = 1.0;

			JPanel menu = new JPanel(new GridBagLayout());
			menu.setBorder(BorderFactory.createEmptyBorder());

			JScrollPane menuScrollPane = new JScrollPane(menu);


			menuScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
			menuScrollPane.setBorder(BorderFactory.createEmptyBorder());

			GridBagConstraints c = new GridBagConstraints();
			c.weighty = 0.0;
			c.weightx = 1.0;
			c.anchor = GridBagConstraints.PAGE_START;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.gridx = 0;
			c.gridy = 0;

			setTrackSwitchesEnabled(false);
			for (JCheckBox trackSwitch : trackSwitches.keySet()) {
				trackSwitch.addActionListener(this);
				menu.add(trackSwitch, c);
				c.gridy++;
			}

			// coverage type
			coverageTypeLabel.setEnabled(false);
			c.insets.set(10,0,0,0);
			menu.add(coverageTypeLabel, c);
			c.gridy++;
			c.insets.set(0,0,0,0);
			coverageTypeBox = new JComboBox(new String[] {COVERAGE_NONE, COVERAGE_TOTAL, COVERAGE_STRAND});
			coverageTypeBox.setSelectedItem(COVERAGE_TOTAL);
			coverageTypeBox.setEnabled(false);
			coverageTypeBox.addActionListener(this);
			menu.add(coverageTypeBox, c);

			// coverage scale
			c.gridy++;
			coverageScaleLabel.setEnabled(false);
			c.insets.set(10,0,0,0);
			menu.add(coverageScaleLabel, c);
			c.gridy++;
			c.insets.set(0,0,0,0);
			coverageScaleBox = new JComboBox(GenomePlot.ReadScale.values());
			coverageScaleBox.setEnabled(false);
			coverageScaleBox.addActionListener(this);
			menu.add(coverageScaleBox, c);

			c.gridy++;
			showFullHeightBox = new JCheckBox("Show all reads", false);
			showFullHeightBox.setEnabled(false);
			showFullHeightBox.addActionListener(this);
			menu.add(showFullHeightBox, c);

			optionsPanel.add(menu, oc);
		}
		return optionsPanel;
	}

	private JPanel getExternalLinkPanel() {
		if (linksPanel == null) { 
			linksPanel = new JPanel(new GridBagLayout());
			linksPanel.setBorder(VisualConstants.createSettingsPanelSubPanelBorder("External links"));

			ensemblLink = LinkUtil.createLink("Ensembl", new AbstractAction() {
				@Override
				public void actionPerformed(ActionEvent arg0) {			
					openExternalBrowser(AnnotationType.ENSEMBL_BROWSER_URL);
				}
			});


			ucscLink = LinkUtil.createLink("UCSC", new AbstractAction() {
				@Override
				public void actionPerformed(ActionEvent arg0) {			
					openExternalBrowser(AnnotationType.UCSC_BROWSER_URL);
				}
			});

			ensemblLink.setEnabled(false);
			ucscLink.setEnabled(false);

			GridBagConstraints c = new GridBagConstraints();
			c.gridx = 0;
			c.gridy = 0;
			//c.gridwidth = 3;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.anchor = GridBagConstraints.NORTHWEST;
			c.weightx = 0;
			c.weighty = 0;

			List<JXHyperlink> importLinks = new LinkedList<JXHyperlink>();
			importLinks.add(ensemblLink);
			importLinks.add(ucscLink);

			final int MAX_ROW_CHARS = 33;

			LinkUtil.addLinks("View this region in *** or *** genome browser.", 
					importLinks, null, c, linksPanel, MAX_ROW_CHARS, null);

			// Panels to take rest of space
			JPanel bottomPanel = new JPanel();
			JPanel rightPanel = new JPanel();

			c.weightx = 0.0;
			c.weighty = 1.0;
			c.fill = GridBagConstraints.VERTICAL;
			c.gridx = 1;
			c.gridy++;
			linksPanel.add(bottomPanel, c);
			c.weightx = 1.0;
			c.weighty = 0.0;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.gridx = 2;
			c.gridy = 1;
			linksPanel.add(rightPanel, c);

			// linksPanel.setMinimumSize(new Dimension(0, 0));
			//			this.setPreferredSize(new Dimension(VisualConstants.LEFT_PANEL_WIDTH, VisualConstants.TREE_PANEL_HEIGHT));

		}

		return linksPanel;
	}


	private void setExternalLinksEnabled() {

		boolean hasLocation = plot != null && plot.getDataView() != null && plot.getDataView().getBpRegion() != null;

		ensemblLink.setEnabled(hasLocation && getExternalLinkUrl(AnnotationType.ENSEMBL_BROWSER_URL).length() > 0);
		ucscLink.setEnabled(hasLocation && getExternalLinkUrl(AnnotationType.UCSC_BROWSER_URL).length() > 0);
	}

	private String getExternalLinkUrl(AnnotationType browser) {
		Genome genome = (Genome) genomeBox.getSelectedItem();
		URL url = annotationManager.getAnnotation(genome, browser).getUrl();

		if (url != null) {
			return url.toString();
		} else {
			return "";
		}

	}

	public void openExternalBrowser(AnnotationType browser) {

		String url = getExternalLinkUrl(browser);	
		Region region = this.plot.getDataView().getBpRegion();
		url = url.replace(AnnotationManager.CHR_LOCATION, region.start.chr.toNormalisedString());
		url = url.replace(AnnotationManager.START_LOCATION, region.start.bp.toString());
		url = url.replace(AnnotationManager.END_LOCATION, region.end.bp.toString());

		try {
			BrowserLauncher.openURL(url);
		} catch (Exception e) {
			application.reportException(e);
		}
	}

	public JPanel createSettingsPanel() {

		settingsPanel.setLayout(new GridBagLayout());

		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 0;
		c.gridx = 0;
		c.insets.set(0, 5, 15, 5);
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weighty = 0;
		c.weightx = 1.0;

		// genome
		settingsPanel.add(getGenomePanel(), c);
		c.gridy++;

		// location
		settingsPanel.add(getLocationPanel(), c);
		c.gridy++;


		// options
		settingsPanel.add(getOptionsPanel(), c);
		c.gridy++;

		// datasets
		c.insets.set(0, 5, 5, 5);
		settingsPanel.add(getDatasetsPanel(), c);
		c.gridy++;

		// external links
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1;
		settingsPanel.add(getExternalLinkPanel(), c);

		return settingsPanel;
	}


	private JPanel getGenomePanel() {
		if (this.genomePanel == null) {
			this.genomePanel = new JPanel(new GridBagLayout());
			genomePanel.setBorder(VisualConstants.createSettingsPanelSubPanelBorder("Genome"));

			GridBagConstraints c = new GridBagConstraints();
			c.gridy = 0;
			c.gridx = 0;
			c.insets.set(5, 0, 5, 0);
			c.anchor = GridBagConstraints.NORTHWEST;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.weighty = 0;
			c.weightx = 1.0;
			c.gridx = 0;

			// genome
			Collection<Genome> genomes = annotationManager.getGenomes();
			for (Genome genome : genomes) {
				genomeBox.addItem(genome);
			}

			// no selection at startup
			genomeBox.setSelectedItem(null);
			genomeBox.addActionListener(this);

			genomePanel.add(genomeBox, c);
		}

		return genomePanel;
	}

	private JPanel getLocationPanel() {
		if (this.locationPanel == null) {

			locationPanel = new JPanel(new GridBagLayout());
			locationPanel.setBorder(VisualConstants.createSettingsPanelSubPanelBorder("Location"));

			GridBagConstraints c = new GridBagConstraints();
			c.gridy = 0;
			c.gridx = 0;
			c.anchor = GridBagConstraints.NORTHWEST;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.weighty = 0;
			c.weightx = 1.0;

			// chromosome
			chrLabel.setEnabled(false);
			locationPanel.add(chrLabel, c);
			c.gridy++;
			chrBox.setEnabled(false);
			c.insets.set(0, 0, 10, 0);
			locationPanel.add(chrBox, c);

			// location
			c.gridy++;
			c.insets.set(0, 0, 0, 0);
			locationLabel.setEnabled(false);
			locationPanel.add(locationLabel, c);
			locationField.setEnabled(false);
			locationField.addActionListener(this);
			c.gridy++;
			c.insets.set(0, 0, 10, 0);
			locationPanel.add(locationField, c);

			// viewsize
			c.gridx = 0;
			c.gridwidth = 5;
			c.gridy++;
			c.insets.set(0, 0, 0, 0);
			viewsizeLabel.setEnabled(false);
			locationPanel.add(viewsizeLabel, c);
			c.gridwidth = 4;
			c.gridy++;
			c.insets.set(0, 0, 10, 0);
			viewsizeField.setEnabled(false);
			viewsizeField.setEditable(false); // view size is never editable
			locationPanel.add(this.viewsizeField, c);

			// go button
			c.gridy++;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.anchor = GridBagConstraints.CENTER;
			goButton.setEnabled(false);
			goButton.addActionListener(this);
			locationPanel.add(goButton, c);
		}
		return locationPanel;
	}

	private void fillChromosomeBox() throws IOException {

		// Gather all chromosome names from all indexed datasets (SAM/BAM)
		TreeSet<String> chromosomeNames = new TreeSet<String>(); 
		for (Interpretation interpretation : interpretations) {
			if (interpretation.type == TrackType.READS) {
				DataBean data = interpretation.primaryData;
				InputStream in = null;
				try {
					in  = data.getContentByteStream();
					chromosomeNames.addAll(SamBamUtils.readChromosomeNames(in));
				} finally {
					IOUtils.closeIfPossible(in);
				}
			}
		}

		// If we still don't have names, go through non-indexed datasets
		if (chromosomeNames.isEmpty()) {
			for (Interpretation interpretation : interpretations) {
				if (interpretation.type == TrackType.REGIONS) {
					DataBean data = interpretation.primaryData;
					File file = Session.getSession().getDataManager().getLocalFile(data);
					List<RegionContent> rows = null;
					try {
						rows = new RegionOperations().loadFile(file);
						for (RegionContent row : rows) {
							chromosomeNames.add(row.region.start.chr.toNormalisedString());
						}
					} catch (URISyntaxException e) {
						e.printStackTrace();
					}
				}
			}
		}

		// Sort them
		LinkedList<Chromosome> chromosomes = new LinkedList<Chromosome>();
		for (String chromosomeName : chromosomeNames) {
			chromosomes.add(new Chromosome(chromosomeName));
		}
		Collections.sort(chromosomes);

		// Fill in the box
		for (Chromosome chromosome : chromosomes) {
			chrBox.addItem(chromosome);
		}
	}

	protected JComponent getColorLabel() {
		return new JLabel("Color: ");
	}

	/**
	 * Changes the visibility of some tracks. This is useful when track is hidden to free some 
	 * screen estate or drawing performance, but the data is still kept in memory (because it's 
	 * needed elsewhere or we just don't care about it).
	 * 
	 * Currently this is used to change between different views of user's bam-files.
	 * 
	 * For more radical changes use updateTracks(), which removes the dataLayers as well.
	 * 
	 */
	public void updateVisibilityForTracks() {
		for (Track track : tracks) {
			if (track.trackGroup != null) {
				for (JCheckBox trackSwitch : trackSwitches.keySet()) {
					track.trackGroup.showOrHide(trackSwitches.get(trackSwitch), trackSwitch.isSelected());
				}

				if (track.trackGroup instanceof ReadTrackGroup) {
					if (coverageTypeBox.getSelectedItem().equals(COVERAGE_NONE)) {
						track.trackGroup.showOrHide("ProfileSNPTrack", false);
						track.trackGroup.showOrHide("ProfileTrack", false);
						track.trackGroup.showOrHide("ReadOverview", false);
					} else 	if (coverageTypeBox.getSelectedItem().equals(COVERAGE_TOTAL)) {
						track.trackGroup.showOrHide("ProfileSNPTrack", true);	
						track.trackGroup.showOrHide("ProfileTrack", false);
						track.trackGroup.showOrHide("ReadOverview", true);
					} else 	if (coverageTypeBox.getSelectedItem().equals(COVERAGE_STRAND)) {
						track.trackGroup.showOrHide("ProfileSNPTrack", false);
						track.trackGroup.showOrHide("ProfileTrack", true);
						track.trackGroup.showOrHide("ReadOverview", true);
					}
				}
			}
		}
	}

	/**
	 * A method defined by the ActionListener interface. Allows this panel to
	 * listen to actions on its components.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == goButton || source == locationField) {

			// disable changing of the genome
			this.genomeBox.setEnabled(false);
			if (!initialised) {

				application.runBlockingTask("initialising genome browser", new Runnable() {
					@Override
					public void run() {
						try {
							// Preload datas in background thread
							initialiseUserDatas();

							// Update UI in Event Dispatch Thread
							SwingUtilities.invokeAndWait(new Runnable() {
								@Override
								public void run() {
									try {
										// Show visualisation
										showVisualisation();
									} catch (Exception e) {
										application.reportException(e);
									}
								}
							});
							
							// Update UI in Event Dispatch Thread, update location only after this block task 
							// has quit to be able to show gene search blocking task. It isn't critical if the timing
							// fails, search will still work, but there the fancy white glass pane won't show.
							SwingUtilities.invokeLater(new Runnable() {
								@Override
								public void run() {

									try {
										updateLocation();
										setExternalLinksEnabled();

										initialised = true;

									} catch (Exception e) {
										application.reportException(e);
									}
								}
							});

						} catch (Exception e) {
							throw new RuntimeException(e);
						}
					}
				});

			} else {

				// Move to correct location
				updateLocation();

				this.setExternalLinksEnabled();
			}

		} else if (datasetSwitches.contains(source) && this.initialised) {

			showFullHeightBox.setSelected(false);
			setFullHeight(false);

			updateTracks();

		} else if (source == coverageScaleBox && this.initialised) {

			updateCoverageScale();

		} else if ((trackSwitches.keySet().contains(source) || source == coverageTypeBox) && this.initialised) {
			updateVisibilityForTracks();
		} 

		// genome selected
		else if (source == genomeBox) {

			Genome genome = (Genome) genomeBox.getSelectedItem();

			// dialog for downloading annotations if not already local
			if (!annotationManager.hasLocalAnnotations(genome)) {
				annotationManager.openDownloadAnnotationsDialog(genome);
			}

			// enable other settings
			this.goButton.setEnabled(true);
			this.chrLabel.setEnabled(true);
			this.chrBox.setEnabled(true);
			this.locationLabel.setEnabled(true);
			this.locationField.setEnabled(true);
			this.viewsizeLabel.setEnabled(true);
			this.viewsizeField.setEnabled(true);
			this.showFullHeightBox.setEnabled(true);

			for (Track track : tracks) {
				track.checkBox.setEnabled(true);
			}

			coverageTypeLabel.setEnabled(true);
			coverageTypeBox.setEnabled(true);

			coverageScaleLabel.setEnabled(true);
			coverageScaleBox.setEnabled(true);

			this.setTrackSwitchesEnabled(true);

		} else if (source == showFullHeightBox && this.initialised) {

			setFullHeight(showFullHeightBox.isSelected());
		}
	}

	private void setFullHeight(boolean fullHeight) {

		if (fullHeight) {
			verticalScroller.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		} else {
			verticalScroller.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_NEVER);
		}

		plot.setFullHeight(fullHeight);		
	}

	private void setTrackSwitchesEnabled(boolean enabled) {
		for (JCheckBox trackSwitch : trackSwitches.keySet()) {
			trackSwitch.setEnabled(enabled);
		}
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		return getVisualisation(Arrays.asList(new DataBean[] { data }));
	}


	@Override
	public JComponent getVisualisation(java.util.List<DataBean> datas) throws Exception {
		


		this.interpretations = interpretUserDatas(datas);

		// List available chromosomes from user data files
		fillChromosomeBox();

		// We can create tracks now that we know the data
		this.tracks.clear();
		createAvailableTracks(); 

		// Create panel with card layout and put message panel there
		JPanel waitPanel = new JPanel(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();

		waitPanel.add(new JLabel("<html><p style=\"font-size: larger\">Please select genome and click " + goButton.getText() + "</p></html>"), c);
		plotPanel.add(waitPanel, WAITPANEL);

		return plotPanel;
	}

	private void updateCoverageScale() {
		// Set scale of profile track containing reads information
		this.plot.setReadScale((ReadScale) this.coverageScaleBox.getSelectedItem());
	}

	private Genome getGenome() {
		return (Genome) genomeBox.getSelectedItem();
	}

	/**
	 * Removes all tracks and data layers and creates new tracks according to current settings.
	 * This is useful when the dataset selection is changed and only those datasets are kept 
	 * in memory that are currently in use.
	 * 
	 * @See updateVisibilityForTracks()
	 */
	private void updateTracks() {

		plot.getDataView().clean();

		Genome genome = getGenome();

		// Add selected annotation tracks
		for (Track track : tracks) {
			if (track.checkBox.isSelected()) {
				switch (track.interpretation.type) {
				case CYTOBANDS:

					URL cytobandUrl = annotationManager.getAnnotation(
							genome, AnnotationManager.AnnotationType.CYTOBANDS).getUrl();

					try {
						CytobandDataSource cytobandDataSource;
						cytobandDataSource = new CytobandDataSource(cytobandUrl);

						TrackFactory.addCytobandTracks(plot, cytobandDataSource);

						this.viewLimiter = new ViewLimiter(plot.getOverviewView().getQueueManager(), 
								cytobandDataSource, plot.getOverviewView());
						this.plot.getDataView().setViewLimiter(viewLimiter);
						this.plot.getOverviewView().setViewLimiter(viewLimiter);

					} catch (FileNotFoundException e) {
						application.reportException(e);
					} catch (URISyntaxException e) {
						application.reportException(e);
					}

					break;

				case GENES:
					// Start 3D effect
					plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, true));

					URL gtfUrl = annotationManager.getAnnotation(
							genome, AnnotationManager.AnnotationType.GTF_TABIX).getUrl();

					URL gtfIndexUrl = annotationManager.getAnnotation(
							genome, AnnotationManager.AnnotationType.GTF_TABIX_INDEX).getUrl();

					URL repeatUrl = annotationManager.getAnnotation(
							genome, AnnotationManager.AnnotationType.REPEAT).getUrl();

					URL repeatIndexUrl = annotationManager.getAnnotation(
							genome, AnnotationManager.AnnotationType.REPEAT_INDEX).getUrl();

					TabixDataSource gtfDataSource;

					try {
						gtfDataSource = new TabixDataSource(gtfUrl, gtfIndexUrl, 
								GtfTabixHandlerThread.class);

						TabixDataSource repeatDataSource = new TabixDataSource(repeatUrl, repeatIndexUrl, BedTabixHandlerThread.class);

						TrackGroup geneGroup = TrackFactory.addGeneTracks(plot, gtfDataSource, repeatDataSource);

						track.setTrackGroup(geneGroup);

					} catch (URISyntaxException e) {
						application.reportException(e);
					} catch (IOException e) {
						application.reportException(e);
					}
					break;

				case REFERENCE:
					// integrated into peaks
					break;

				case TRANSCRIPTS:
					// integrated into genes
					break;
				}
			}
		}

		// Add selected read tracks
		for (Track track : tracks) {
			if (track.checkBox.isSelected()) {

				File file;
				try {
					file = track.interpretation.primaryData == null ? null : Session
							.getSession().getDataManager().getLocalFile(
									track.interpretation.primaryData);
					DataSource treatmentData;
					if (track.interpretation.type == TrackType.READS) {

						URL fastaUrl = annotationManager.getAnnotation(
								genome, AnnotationManager.AnnotationType.REFERENCE).getUrl();

						URL fastaIndexUrl = annotationManager.getAnnotation(
								genome, AnnotationManager.AnnotationType.REFERENCE_INDEX).getUrl();

						IndexedFastaDataSource refSeqDataSource = new IndexedFastaDataSource(fastaUrl, fastaIndexUrl);

						if (track.interpretation.summaryDatas.size() == 0) {
							// No precomputed summary data
							TrackFactory.addThickSeparatorTrack(plot);
							treatmentData = createReadDataSource(track.interpretation.primaryData, track.interpretation.indexData, tracks);

							TrackGroup readGroup = TrackFactory.addReadTracks(
									plot, treatmentData, 
									refSeqDataSource, 
									track.interpretation.primaryData.getName());

							track.setTrackGroup(readGroup);

						} else { 
							// Has precomputed summary data
							TrackFactory.addThickSeparatorTrack(plot);
							treatmentData = createReadDataSource(track.interpretation.primaryData, track.interpretation.indexData, tracks);
							TrackGroup readGroupWithSummary = TrackFactory.addReadSummaryTracks(
									plot, treatmentData, refSeqDataSource, 
									track.interpretation.primaryData.getName(), new TabixDataSource(file.toURI().toURL(), null, TabixSummaryHandlerThread.class));
							track.setTrackGroup(readGroupWithSummary);
						}
					}
				} catch (IOException e) {
					application.reportException(e);
				} catch (MicroarrayException e) {
					application.reportException(e);
				} catch (URISyntaxException e) {
					application.reportException(e);
				}
			}
		}

		// Add selected peak tracks
		for (Track track : tracks) {
			if (track.checkBox.isSelected()) {

				URL fileUrl = null;

				if (track.interpretation.primaryData != null) {
					File file;
					try {
						file = Session.getSession().getDataManager().getLocalFile(
								track.interpretation.primaryData);
						fileUrl = file.toURI().toURL();

					} catch (IOException e) {
						application.reportException(e);
					}
				}

				DataSource regionData;
				switch (track.interpretation.type) {
				case REGIONS:
					TrackFactory.addThickSeparatorTrack(plot);
					TrackFactory.addTitleTrack(plot, track.interpretation.primaryData.getName());

					try {
						regionData = new ChunkDataSource(fileUrl, new BEDParserWithCoordinateConversion(), ChunkTreeHandlerThread.class);
						((ChunkDataSource)regionData).checkSorting();
						TrackFactory.addPeakTrack(plot, regionData);

					} catch (FileNotFoundException e) {
						application.reportException(e);
					} catch (URISyntaxException e) {
						application.reportException(e);
					} catch (IOException e) {
						application.reportException(e);
					} catch (MicroarrayException e) {
						application.reportException(e);
					} catch (UnsortedDataException e) {
						application.showDialog("Unsorted data", e.getMessage(), null, Severity.WARNING, true);
					}
					break;
				case REGIONS_WITH_HEADER:
					TrackFactory.addThickSeparatorTrack(plot);
					TrackFactory.addTitleTrack(plot, track.interpretation.primaryData.getName());

					try {
						regionData = new ChunkDataSource(fileUrl, new HeaderTsvParser(), ChunkTreeHandlerThread.class);
						TrackFactory.addPeakTrack(plot, regionData);

					} catch (FileNotFoundException e) {
						application.reportException(e);
					} catch (URISyntaxException e) {
						application.reportException(e);
					}
					break;
				case VCF:
					TrackFactory.addThickSeparatorTrack(plot);
					TrackFactory.addTitleTrack(plot, track.interpretation.primaryData.getName());

					try {
						regionData = new ChunkDataSource(fileUrl, new VcfParser(), ChunkTreeHandlerThread.class);
						TrackFactory.addPeakTrack(plot, regionData);

					} catch (FileNotFoundException e) {
						application.reportException(e);
					} catch (URISyntaxException e) {
						application.reportException(e);
					}
					break;

				}
			}
		}

		// End 3D effect
		plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, false));

		// Set track visibility
		updateVisibilityForTracks();
	}

	private void showVisualisation() {

		//Clean old data layers
		if (plot != null) {
			plot.clean();
		}

		// Create the chart panel with tooltip support				
		TooltipAugmentedChartPanel chartPanel = new TooltipAugmentedChartPanel();
		this.plot = new GenomePlot(chartPanel, true);
		((NonScalableChartPanel)chartPanel).setGenomePlot(plot);

		//Set location to plot to avoid trouble in track initialization. 
		//Can't do this with updateLocation, because it would lose gene search when the 
		//tracks clear all data layers
		plot.getDataView().setBpRegion(new RegionDouble(
				DEFAULT_LOCATION - DEFAULT_VIEWSIZE / 2.0, DEFAULT_LOCATION + DEFAULT_VIEWSIZE / 2.0, 
				(Chromosome)chrBox.getSelectedItem()), true);
		
		updateCoverageScale();
		
		updateTracks();

		// Initialise the plot
		plot.addDataRegionListener(this);

		// Wrap GenomePlot in a panel
		chartPanel.setChart(new JFreeChart(plot));
		chartPanel.setCursor(new Cursor(Cursor.HAND_CURSOR));

		// Add mouse listeners
		for (View view : plot.getViews()) {
			chartPanel.addMouseListener(view);
			chartPanel.addMouseMotionListener(view);
			chartPanel.addMouseWheelListener(view);
		}

		// Add selection listener (but try to remove first old one that would prevent removal of the visualization) 
		application.removeClientEventListener(this);
		application.addClientEventListener(this);

		// Put panel on top of card layout
		if (plotPanel.getComponentCount() == 2) {
			plotPanel.remove(1);
		}

		verticalScroller = new JScrollPane(chartPanel);

		setFullHeight(showFullHeightBox.isSelected());

		plotPanel.add(verticalScroller, PLOTPANEL);
		plotPanel.addComponentListener(this);
		CardLayout cl = (CardLayout) (plotPanel.getLayout());
		cl.show(plotPanel, PLOTPANEL);
	}

	private GeneIndexActions getGeneIndexActions() {

		if (gia == null) {
			Genome genome = getGenome();

			// Create gene name index
			try {

				URL gtfUrl = annotationManager.getAnnotation(
						genome, AnnotationManager.AnnotationType.GTF_TABIX).getUrl();

				URL gtfIndexUrl = annotationManager.getAnnotation(
						genome, AnnotationManager.AnnotationType.GTF_TABIX_INDEX).getUrl();

				URL geneUrl = annotationManager.getAnnotation(
						genome, AnnotationManager.AnnotationType.GENE_CHRS).getUrl();


				TabixDataSource gtfDataSource = new TabixDataSource(gtfUrl, gtfIndexUrl, GtfTabixHandlerThread.class);
				LineDataSource geneDataSource = new LineDataSource(geneUrl, GeneSearchHandler.class);

				gia = new GeneIndexActions(plot.getDataView().getQueueManager(), gtfDataSource, geneDataSource);

			} catch (Exception e) {
				application.reportException(e);
			}
		}
		return gia;
	}

	private void initialiseUserDatas() throws IOException {
		for (Interpretation interpretation : interpretations) {
			initialiseUserData(interpretation.primaryData);
			initialiseUserData(interpretation.indexData);
			for (DataBean summaryData : interpretation.summaryDatas) {
				initialiseUserData(summaryData);
			}
		}
	}

	private void initialiseUserData(DataBean data) throws IOException {
		if (data != null) {
			Session.getSession().getDataManager().getLocalFile(data);
		}
	}


	/**
	 * Create DataSource either for SAM/BAM or ELAND data files.
	 * 
	 * @param tracks
	 * 
	 * @param file
	 * @return
	 * @throws MicroarrayException
	 *             if index file is not selected properly
	 * @throws IOException
	 *             if opening data files fails
	 * @throws URISyntaxException 
	 */
	public DataSource createReadDataSource(DataBean data, DataBean indexData, List<Track> tracks)
			throws MicroarrayException, IOException, URISyntaxException {
		DataSource dataSource = null;

		// Convert data bean into file
		File file = data == null ? null : Session.getSession().getDataManager().getLocalFile(data);

		URL fileUrl = file.toURI().toURL();

		if (data.getName().contains(".bam-summary")) {
			dataSource = new TabixSummaryDataSource(fileUrl);

		} else if (data.getName().contains(".bam") || data.getName().contains(".sam")) {
			File indexFile = Session.getSession().getDataManager().getLocalFile(indexData);
			URL indexFileUrl = indexFile.toURI().toURL();
			dataSource = new SAMDataSource(fileUrl, indexFileUrl);

		} else {
			dataSource = new ChunkDataSource(fileUrl, new ElandParser(), ChunkTreeHandlerThread.class);
		}

		return dataSource;
	}

	private boolean isIndexData(DataBean bean) {
		return bean.getName().endsWith(".bai");
	}

	@Override
	public boolean canVisualise(DataBean data) throws MicroarrayException {
		return canVisualise(Arrays.asList(new DataBean[] { data }));
	}

	@Override
	public boolean canVisualise(java.util.List<DataBean> datas) throws MicroarrayException {
		return interpretUserDatas(datas) != null;
	}

	public class ObjVariable extends Variable {

		public Object obj;

		public ObjVariable(Object obj) {
			super(null, null);
			this.obj = obj;
		}
	}

	public void regionChanged(Region bpRegion) {
		setCoordinateFields(bpRegion.getMid(), bpRegion.getLength());
		this.lastLocation = bpRegion.getMid();
		this.lastViewsize = bpRegion.getLength();
	}

	private List<Interpretation> interpretUserDatas(List<DataBean> datas) {
		LinkedList<Interpretation> interpretations = new LinkedList<Interpretation>();

		// Find interpretations for all primary data types
		for (DataBean data : datas) {
			


			if (data.isContentTypeCompatitible("text/plain")) {
				// ELAND result / export
				interpretations.add(new Interpretation(TrackType.READS, data));

			} else if (data.isContentTypeCompatitible("text/bed")) {
				// BED (ChIP-seq peaks)
				interpretations.add(new Interpretation(TrackType.REGIONS, data));

			} else if (data.isContentTypeCompatitible("text/tab")) {
				// peaks (with header in the file)
				interpretations.add(new Interpretation(TrackType.REGIONS_WITH_HEADER, data));

			} else if ((data.isContentTypeCompatitible("application/bam"))) {
				// BAM file
				interpretations.add(new Interpretation(TrackType.READS, data));

			} else if ((data.isContentTypeCompatitible("text/vcf"))) {
				// Vcf file
				interpretations.add(new Interpretation(TrackType.VCF, data));
			}
		}

		// Find interpretations for all secondary data types
		for (DataBean data : datas) {

			// Find the interpretation to add this secondary data to
			Interpretation primaryInterpretation = null;
			for (Interpretation interpretation : interpretations) {
				if (data.getName().startsWith(interpretation.primaryData.getName())) {
					primaryInterpretation = interpretation;
					break;
				}
			}

			if (primaryInterpretation == null) {
				return null; // could not bound this secondary data to any primary data
			}

			if ((data.isContentTypeCompatitible("application/octet-stream")) &&
					(data.getName().contains(".bam-summary"))) {
				// BAM summary file (from custom preprocessor)
				primaryInterpretation.summaryDatas.add(data);

			} else if ((data.isContentTypeCompatitible("application/octet-stream")) &&
					(isIndexData(data))) {
				// BAI file
				if (primaryInterpretation.indexData != null) {
					return null; // already taken, could not bind this secondary data to any primary data
				}
				primaryInterpretation.indexData = data;
			}
		}

		// Check that interpretations are now fully initialised
		for (Interpretation interpretation : interpretations) {
			if (interpretation.primaryData.getName().endsWith(".bam") && interpretation.indexData == null) {
				return null; // BAM is missing BAI
			}
		}


		return interpretations;
	}

	@Override
	public boolean isForSingleData() {
		return true;
	}

	@Override
	public boolean isForMultipleDatas() {
		return true;
	}

	/**
	 * Update genome browser to location given in the location panel.
	 * 
	 */
	private void updateLocation() {

		boolean isSearch = false;

		if (!locationField.getText().isEmpty() && !GeneIndexActions.checkIfNumber(locationField.getText())) {
			// If gene name was given, search for it

			requestGeneSearch();
			isSearch = true;
		}

		// Fill in initial position if not filled in
		if (locationField.getText().trim().isEmpty() || isSearch) {

			setCoordinateFields(DEFAULT_LOCATION, null);
		}

		if (viewsizeField.getText().trim().isEmpty()) {

			setCoordinateFields(null, DEFAULT_VIEWSIZE);
			lastViewsize = DEFAULT_VIEWSIZE;
		}


		plot.moveDataBpRegion((Chromosome) chrBox.getSelectedItem(),
				Long.parseLong(locationField.getText()), lastViewsize);

		// Set scale of profile track containing reads information
		this.plot.setReadScale((ReadScale) this.coverageScaleBox.getSelectedItem());
	}

	/**
	 * Null keeps the existing content.
	 * 
	 * @param location
	 * @param viewsize
	 */
	private void setCoordinateFields(Long location, Long viewsize) {
		if (location != null) {
			locationField.setText(location.toString());
		}

		if (viewsize != null) {
			if (viewsize > 1000000) {
				viewsizeField.setText(Math.round(((float)viewsize) / 1000000f) + " Mb");
			} else if (viewsize > 1000) {
				viewsizeField.setText(Math.round(((float)viewsize) / 1000f) + " kb");
			} else {
				viewsizeField.setText(viewsize + "");
			}
		}
	}

	private void requestGeneSearch() {

		application.runBlockingTask("searching gene", new Runnable() {

			@Override
			public void run() {

				int TIME_OUT = 30*1000;
				int INTERVAL = 100;

				long startTime = System.currentTimeMillis();

				while (System.currentTimeMillis() < startTime + TIME_OUT) {

					if (!geneSearchDone) {
						try {
							Thread.sleep(INTERVAL);
						} catch (InterruptedException e) {
							//Just continue
						}
					} else {
						break;
					}
				}		

				if (geneSearchDone) {

					geneSearchDone = false;

				} else {

					//Give up
					SwingUtilities.invokeLater(new Runnable() {

						@Override
						public void run() {

							application.showDialog("Search failed",
									"Unexpected error happened in the search. Please inform the developers if the problem persists.", null,
									Severity.WARNING, true,
									DetailsVisibility.DETAILS_ALWAYS_HIDDEN, null);
						}

					});
				}
			}
		});

		getGeneIndexActions().requestLocation(locationField.getText(), new GeneIndexActions.GeneLocationListener() {

			@Override
			public void geneLocation(Region geneLocation) {

				geneSearchDone = true;

				if (geneLocation == null) {

					// Move to last known location
					updateLocation();

					// Tell the user 
					application.showDialog("Not found",
							"Gene was not found", null,
							Severity.INFO, true,
							DetailsVisibility.DETAILS_ALWAYS_HIDDEN, null);

				} else {

					// Update coordinate controls with gene's location

					Chromosome resultChr = new Chromosome(geneLocation.start.chr);

					chrBox.setSelectedItem(resultChr);

					if (chrBox.getSelectedItem().equals(resultChr)) {

						setCoordinateFields((geneLocation.end.bp + geneLocation.start.bp) / 2, (geneLocation.end.bp - geneLocation.start.bp) * 2);
						updateLocation();
					} else {
						application.showDialog("Different chromosome", 
								"Searched gene was found from chromosome " + resultChr + " but there is no data for that chromosome", "" + geneLocation, 
								Severity.INFO, true, 
								DetailsVisibility.DETAILS_HIDDEN, null);
					}
				}
			}
		});
	}

	@Override
	public void componentHidden(ComponentEvent arg0) {
		// Ignore
	}

	@Override
	public void componentMoved(ComponentEvent arg0) {
		// Ignore
	}

	@Override
	public void componentResized(ComponentEvent arg0) {
		this.updateLocation();
		plot.redraw();
	}

	@Override
	public void componentShown(ComponentEvent arg0) {
		// Ignore
	}

	@Override
	public void propertyChange(PropertyChangeEvent event) {
		if (event instanceof PointSelectionEvent) {

			IntegratedEntity sel = application.getSelectionManager().getSelectionManager(null).getPointSelection();

			// Check if we can process this
			if (sel.containsKey("chromosome") && sel.containsKey("start")) {

				// Move to selected region 
				chrBox.setSelectedItem(new Chromosome(sel.get("chromosome")));
				long start = Long.parseLong(sel.get("start"));

				long end = -1;
				if (sel.containsKey("end")) {
					end = Long.parseLong(sel.get("end"));
				} else {
					end = start;
				}
				setCoordinateFields((end + start) / 2, (end - start) * 2);
			}

			// Update
			updateLocation();
		}
	}

	@Override
	public void removeVisualisation() {

		super.removeVisualisation();

		plotPanel.removeComponentListener(this);
		plotPanel.removeAll();

		if (plot != null) {
			plot.clean();
			plot = null;
		}

		//Remove references to tracks and data to free memory, even if the (hidden) parameter panel keeps actionListener
		//references to this object preventing garbage collection (when visualization is changed to none)
		if (tracks != null) {
			tracks.clear();
		}
		gia = null;

		application.removeClientEventListener(this);

	}
}
