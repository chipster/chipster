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

import org.jfree.chart.JFreeChart;

import fi.csc.chipster.tools.gbrowser.SamBamUtils;
import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.selection.PointSelectionEvent;
import fi.csc.microarray.client.visualisation.NonScalableChartPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomePlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GtfHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.HeaderTsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.GenomeAnnotation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.gbrowser.index.GeneIndexActions;
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
		HIDDEN(false);

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

		public TrackGroup getTrackGroup() {
			return trackGroup;
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
	private JPanel plotPanel = new JPanel(new CardLayout());

	private JButton goButton = new JButton("Go");

	private JLabel locationLabel = new JLabel("Location (gene or position)");
	private JTextField locationField = new JTextField();

	private JLabel viewsizeLabel = new JLabel("View size");
	private JTextField viewsizeField = new JTextField();

	private JLabel chrLabel = new JLabel("Chromosome");
	private JComboBox chrBox = new JComboBox();

	private JComboBox genomeBox = new JComboBox();

	private Object visibleChromosome;

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


	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);

		// initialize annotations
		this.annotationManager = new AnnotationManager();
		this.annotationManager.initialize();

		trackSwitches.put(new JCheckBox("Reads", true), "Reads");
		//		trackSwitches.put(new JCheckBox("Highlight SNPs", false), "highlightSNP");
		//		trackSwitches.put(new JCheckBox("Coverage and SNPs", true), "ProfileSNPTrack");
		//		trackSwitches.put(new JCheckBox("Strand-specific coverage", false), "ProfileTrack");

		//		trackSwitches.put(new JCheckBox("Quality coverage", false), "QualityCoverageTrack"); // TODO re-enable quality coverage
		trackSwitches.put(new JCheckBox("Density graph", false), "GelTrack");
		//trackSwitches.put(new JCheckBox("Low complexity regions", false), "RepeatMaskerTrack"); // TODO re-enable dbSNP view
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
		tracks.add(new Track(AnnotationManager.AnnotationType.ANNOTATIONS.getId(), new Interpretation(TrackType.GENES, null)));
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
			showFullHeightBox = new JCheckBox("Show full height", false);
			showFullHeightBox.setEnabled(false);
			showFullHeightBox.addActionListener(this);
			menu.add(showFullHeightBox, c);

			optionsPanel.add(menu, oc);
		}
		return optionsPanel;
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

		c.fill = GridBagConstraints.BOTH;

		// options
		settingsPanel.add(getOptionsPanel(), c);
		c.gridy++;

		// datasets
		c.weighty = 0.5;
		c.insets.set(0, 5, 5, 5);
		settingsPanel.add(getDatasetsPanel(), c);

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
						track.trackGroup.showOrHide("highlightSNP", true);
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

									// Show visualisation
									showVisualisation();
									
									updateLocation();

									// Create tracks only once
									initialised = true;

									// Set track visibility
									updateVisibilityForTracks();
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
			}

		} else if ((datasetSwitches.contains(source) || source == coverageScaleBox) && this.initialised) {

			showFullHeightBox.setSelected(false);
			setFullHeight(false);

			showVisualisation();
			updateVisibilityForTracks();	        	        

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

	private void showVisualisation() {

		//Clean old data layers
		if (plot != null) {
			plot.clean();
		}

		try {

			// Create the chart panel with tooltip support				
			TooltipAugmentedChartPanel chartPanel = new TooltipAugmentedChartPanel();
			this.plot = new GenomePlot(chartPanel, true);
			((NonScalableChartPanel)chartPanel).setGenomePlot(plot);

			// Set scale of profile track containing reads information
			this.plot.setReadScale((ReadScale) this.coverageScaleBox.getSelectedItem());

			Genome genome = (Genome) genomeBox.getSelectedItem();
			initGeneIndex(genome);

			// Add selected annotation tracks
			for (Track track : tracks) {
				if (track.checkBox.isSelected()) {
					switch (track.interpretation.type) {
					case CYTOBANDS:

						URL cytobandUrl = annotationManager.getAnnotation(
								genome, AnnotationManager.AnnotationType.CYTOBANDS).getUrl();
						URL regionsUrl = annotationManager.getAnnotation(
								genome, AnnotationManager.AnnotationType.CYTOBANDS_SEQ_REGION).getUrl();
						URL coordUrl = annotationManager.getAnnotation(
								genome, AnnotationManager.AnnotationType.CYTOBANDS_COORD_SYSTEM).getUrl();

						CytobandDataSource cytobandDataSource = new CytobandDataSource(cytobandUrl, regionsUrl, coordUrl);

						TrackFactory.addCytobandTracks(plot, cytobandDataSource);

						this.viewLimiter = new ViewLimiter(plot.getOverviewView().getQueueManager(), 
								cytobandDataSource, plot.getOverviewView());
						this.plot.getDataView().setViewLimiter(viewLimiter);

						break;

					case GENES:
						// Start 3D effect
						plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, true));

						URL gtfUrl = annotationManager.getAnnotation(
								genome, AnnotationManager.AnnotationType.ANNOTATIONS).getUrl();

						LineDataSource gtfDataSource = new LineDataSource(gtfUrl, GtfHandlerThread.class);

						TrackGroup geneGroup = TrackFactory.addGeneTracks(plot, gtfDataSource);

						track.setTrackGroup(geneGroup);
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

					File file = track.interpretation.primaryData == null ? null : Session
							.getSession().getDataManager().getLocalFile(
									track.interpretation.primaryData);
					DataSource treatmentData;
					if (track.interpretation.type == TrackType.READS) {

						FastaDataSource refSeqDataSource = new FastaDataSource();

						for (GenomeAnnotation annotation : annotationManager.getAnnotations(
								genome, AnnotationManager.AnnotationType.REFERENCE)) {

							refSeqDataSource.put(annotation.chr, annotation.getUrl());
						}

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
									track.interpretation.primaryData.getName(), new TabixDataSource(file.toURI().toURL()));
							track.setTrackGroup(readGroupWithSummary);
						}
					}
				}
			}

			// Add selected peak tracks
			for (Track track : tracks) {
				if (track.checkBox.isSelected()) {

					URL fileUrl = null;

					if (track.interpretation.primaryData != null) {
						File file = Session.getSession().getDataManager().getLocalFile(
								track.interpretation.primaryData);
						fileUrl = file.toURI().toURL();
					}

					DataSource peakData;
					switch (track.interpretation.type) {
					case REGIONS:
						TrackFactory.addThickSeparatorTrack(plot);
						TrackFactory.addTitleTrack(plot, track.interpretation.primaryData.getName());
						peakData = new ChunkDataSource(fileUrl, new BEDParser(), ChunkTreeHandlerThread.class);
						TrackFactory.addPeakTrack(plot, peakData);
						break;
					case REGIONS_WITH_HEADER:
						TrackFactory.addThickSeparatorTrack(plot);
						TrackFactory.addTitleTrack(plot, track.interpretation.primaryData.getName());
						peakData = new ChunkDataSource(fileUrl, new HeaderTsvParser(), ChunkTreeHandlerThread.class);
						TrackFactory.addHeaderPeakTrack(plot, peakData);
						break;
					}
				}
			}

			// End 3D effect
			plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, false));

			//resetLocationFields();

			// Initialise the plot
			plot.addDataRegionListener(this);

			//updateLocation();
			move();

			// Remember chromosome
			visibleChromosome = chrBox.getSelectedItem();

			// Wrap GenomePlot in a panel
			chartPanel.setChart(new JFreeChart(plot));
			chartPanel.setCursor(new Cursor(Cursor.HAND_CURSOR));

			// Add mouse listeners
			for (View view : plot.getViews()) {
				chartPanel.addMouseListener(view);
				chartPanel.addMouseMotionListener(view);
				chartPanel.addMouseWheelListener(view);
			}

			// Add selection listener
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

		} catch (Exception e) {
			application.reportException(e);
		}
	}

	private void resetLocationFieldsIfEmpty() {
		
		// Fill in initial position if not filled in
		if (locationField.getText().trim().isEmpty()) {
			
			updateCoordinateFields(DEFAULT_LOCATION, null);
		}
		
		if (viewsizeField.getText().trim().isEmpty()) {
			
			updateCoordinateFields(null, DEFAULT_VIEWSIZE);
			lastViewsize = DEFAULT_VIEWSIZE;
		}
	}

	private void initGeneIndex(Genome genome) {

		// Create gene name index
		gia = null;
		try {

			URL gtfUrl = annotationManager.getAnnotation(
					genome, AnnotationManager.AnnotationType.ANNOTATIONS).getUrl();

			LineDataSource gtfDataSource = new LineDataSource(gtfUrl, GtfHandlerThread.class);

			gia = new GeneIndexActions(plot.getDataView().getQueueManager(), gtfDataSource);

		} catch (Exception e) {
			application.reportException(e);
		}
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

		if (file.getName().contains(".bam-summary")) {
			dataSource = new TabixDataSource(fileUrl);

		} else if (file.getName().contains(".bam") || file.getName().contains(".sam")) {
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
		updateCoordinateFields(bpRegion.getMid(), bpRegion.getLength());
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

			} else if ((data.isContentTypeCompatitible("application/octet-stream")) &&
					(data.getName().endsWith(".bam"))) {
				// BAM file
				interpretations.add(new Interpretation(TrackType.READS, data));
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
	 * If chromosome changes, reinitialize everything (because some old
	 * information is left inside the tracks). Otherwise, simply move currently
	 * viewed bp region.
	 * 
	 * TODO Instead of showVisualisation, clean track contents. This is nicer
	 * because we don't have to reinitialize the tracks and track group options
	 * are saved.
	 */
	private void updateLocation() {

		try {

			if (!GeneIndexActions.checkIfNumber(locationField.getText())) {
				// If gene name was given, search for it
				if (!locationField.getText().equals("")) {
					requestGeneSearch();
				}
				locationField.setText("");
				move(); //Init view.bpRegion if necesssary
				
			} else {

				// Check how large the update in location was 
				if (visibleChromosome != null && visibleChromosome != chrBox.getSelectedItem()) {

					// Chromosome changed, redraw everything
					showVisualisation();
					updateVisibilityForTracks();

				} else {
					// Only bp position within chromosome changed, move there
					move();
				}		
			}
		} catch (Exception e) {
			application.reportException(e);
		}
	}

	private void requestGeneSearch() {

		Chromosome chr = (Chromosome) chrBox.getSelectedItem();

		gia.requestLocation(locationField.getText(), chr, new GeneIndexActions.GeneLocationListener() {

			@Override
			public void geneLocation(Region geneLocation) {

				if (geneLocation == null) {

					// Move to last known location
					if (lastLocation != null && lastViewsize != null) {
						updateCoordinateFields(lastLocation, lastViewsize);
					} else {
						updateCoordinateFields(DEFAULT_LOCATION, DEFAULT_VIEWSIZE);
					}

					// Tell the user 
					application.showDialog("Not found",
							"Gene was not found", null,
							Severity.INFO, true,
							DetailsVisibility.DETAILS_ALWAYS_HIDDEN, null);

				} else {

					// Update coordinate controls with gene's location
					chrBox.setSelectedItem(new Chromosome(geneLocation.start.chr));
					updateCoordinateFields((geneLocation.end.bp + geneLocation.start.bp) / 2, (geneLocation.end.bp - geneLocation.start.bp) * 2);
					updateLocation();
				}
			}
		});
	}

	/**
	 * Null keeps the existing content.
	 * 
	 * @param location
	 * @param viewsize
	 */
	private void updateCoordinateFields(Long location, Long viewsize) {
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


	private void move() {
		
		resetLocationFieldsIfEmpty();

		plot.moveDataBpRegion((Chromosome) chrBox.getSelectedItem(),
				Long.parseLong(locationField.getText()), lastViewsize);

		// Set scale of profile track containing reads information
		this.plot.setReadScale((ReadScale) this.coverageScaleBox.getSelectedItem());
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
			if (sel.containsKey("chromosome") && sel.containsKey("start") && sel.containsKey("end")) {

				// Move to selected region 
				chrBox.setSelectedItem(new Chromosome(sel.get("chromosome")));
				long start = Long.parseLong(sel.get("start"));
				long end = Long.parseLong(sel.get("end"));
				updateCoordinateFields((end + start) / 2, (end - start) * 2);
			}

			// Update
			updateLocation();
		}
	}

	@Override
	public void removeVisualisation() {

		super.removeVisualisation();

		if (plot != null) {
			plot.clean();
		}
		// Keeping database consumes maybe 30 MB of RAM per genome, but subsequent GB visualisations start faster
		//		gia.clean();	
		//		this.gia = null;

		application.removeClientEventListener(this);

	}
}
