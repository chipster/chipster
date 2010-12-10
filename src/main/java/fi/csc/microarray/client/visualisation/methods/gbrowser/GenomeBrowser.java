package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.CardLayout;
import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
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

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import fi.csc.chipster.tools.gbrowser.SamBamUtils;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.visualisation.NonScalableChartPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomePlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.SAMHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.CytobandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.GeneParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.HeaderTsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SNPParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SequenceParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TranscriptParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.GenomeAnnotation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.gbrowser.index.GeneIndexActions;
import fi.csc.microarray.util.IOUtils;

/**
 * Chipster style visualisation for genome browser.
 * 
 * @author Petri Klemel�, Aleksi Kallio
 */
public class GenomeBrowser extends Visualisation implements ActionListener,
		RegionListener, FocusListener, ComponentListener {


	private static final String DEFAULT_ZOOM = "100000";
	private static final String DEFAULT_LOCATION = "1000000";
	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";

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
		PEAKS(true),
		PEAKS_WITH_HEADER(true), 
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

	private JLabel zoomLabel = new JLabel("Zoom");
	private JTextField zoomField = new JTextField(10);
	
	private JLabel chrLabel = new JLabel("Chromosome");
	private JComboBox chrBox = new JComboBox();
	
	private JComboBox genomeBox = new JComboBox();
	
	private Object visibleChromosome;

	private AnnotationManager annotationManager;

	private JLabel coverageScaleLabel = new JLabel("Coverage scale");
	private JComboBox coverageScaleBox = new JComboBox();

	private GeneIndexActions gia;

	private boolean initialised;
	
	private Map<JCheckBox, String> trackSwitches = new LinkedHashMap<JCheckBox, String>();
	private Set<JCheckBox> datasetSwitches = new HashSet<JCheckBox>();
	private Long lastLocation;
	private Long lastZoom;
	
	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);

		// initialize annotations
		this.annotationManager = new AnnotationManager();
		this.annotationManager.initialize();
		
		trackSwitches.put(new JCheckBox("Reads", true), "Reads");
		trackSwitches.put(new JCheckBox("Highlight SNPs", false), "highlightSNP");
		trackSwitches.put(new JCheckBox("Coverage and SNP's", true), "ProfileSNPTrack");
		trackSwitches.put(new JCheckBox("Strand-specific coverage", false), "ProfileTrack");
		trackSwitches.put(new JCheckBox("Quality coverage", false), "QualityCoverageTrack");
		trackSwitches.put(new JCheckBox("Density graph", false), "GelTrack");
//		trackSwitches.put(new JCheckBox("Show common SNP's", false), "changeSNP"); // TODO re-enable dbSNP view
	}

	@Override
	public JPanel getParameterPanel() {

		if (paramPanel == null) {
			paramPanel = new JPanel();
			paramPanel.setLayout(new GridBagLayout());
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

			JPanel settings = this.createSettingsPanel();

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Settings", settings);

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
		tracks.add(new Track(AnnotationManager.AnnotationType.GENES.getId(), new Interpretation(TrackType.GENES, null)));
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
		if (this.datasetsPanel == null) { 
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
			
			JScrollPane scrollPane = new JScrollPane(datasetSwitchesPanel);
			scrollPane.setBorder(BorderFactory.createEmptyBorder());

			scrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);		
			datasetsPanel.add(scrollPane, dc);
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

			// coverage scale
			coverageScaleLabel.setEnabled(false);
			c.insets.set(10,0,0,0);
			menu.add(coverageScaleLabel, c);
			c.gridy++;
			c.insets.set(0,0,0,0);
			coverageScaleBox = new JComboBox(GenomePlot.ReadScale.values());
			coverageScaleBox.setEnabled(false);
			coverageScaleBox.addActionListener(this);
			menu.add(coverageScaleBox, c);

			optionsPanel.add(menu, oc);
		}
		return optionsPanel;
	}
	
	public JPanel createSettingsPanel() {

		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

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
//		c.weighty = 1.0;
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

			// zoom
			c.gridx = 0;
			c.gridwidth = 5;
			c.gridy++;
			c.insets.set(0, 0, 0, 0);
			zoomLabel.setEnabled(false);
			locationPanel.add(zoomLabel, c);
			c.gridwidth = 4;
			c.gridy++;
			c.insets.set(0, 0, 10, 0);
			zoomField.setEnabled(false);
			locationPanel.add(this.zoomField, c);
			this.zoomField.addActionListener(this);

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
		
		// Gather all chromosome names from all read datasets
		TreeSet<String> chromosomeNames = new TreeSet<String>(); 
		for (Interpretation interpretation : interpretations) {
			TrackType trackType = interpretation.type;
			if (trackType == TrackType.READS) {
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
			}
		}
	}

	/**
	 * A method defined by the ActionListener interface. Allows this panel to
	 * listen to actions on its components.
	 */
	public void actionPerformed(ActionEvent e) {
		Object source = e.getSource();

		if (source == goButton || source == locationField || source == zoomField) {

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
	        showVisualisation();
	        updateVisibilityForTracks();

		} else if (trackSwitches.keySet().contains(source) && this.initialised) {
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
			this.zoomLabel.setEnabled(true);
			this.zoomField.setEnabled(true);
			
			for (Track track : tracks) {
				track.checkBox.setEnabled(true);
			}
			
			coverageScaleLabel.setEnabled(true);
			coverageScaleBox.setEnabled(true);
			
			this.setTrackSwitchesEnabled(true);
		}		
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

		// Create tracks only once
		initialised = true;

		try {
			
			Genome genome = (Genome) genomeBox.getSelectedItem();
			
			// Create gene name index
			gia = null;
			try {
				gia = GeneIndexActions.getInstance(genome, createAnnotationDataSource(annotationManager.getAnnotation(genome, AnnotationManager.AnnotationType.GENES).getUrl(),	new GeneParser()));
			} catch (Exception e) {
				application.reportException(e);
			}
			
			// Create the plot
			ChartPanel chartPanel = new NonScalableChartPanel();
			this.plot = new GenomePlot(chartPanel, true);
			
			// Set scale of profile track containing reads information
			this.plot.setReadScale((ReadScale) this.coverageScaleBox.getSelectedItem());


			// Add selected annotation tracks
			for (Track track : tracks) {
				if (track.checkBox.isSelected()) {
					switch (track.interpretation.type) {
					case CYTOBANDS:
						TrackFactory.addCytobandTracks(plot,
								createAnnotationDataSource(
										annotationManager.getAnnotation(
												genome, AnnotationManager.AnnotationType.CYTOBANDS).getUrl(),
										new CytobandParser()));
						break;
						
					case GENES:
						// Start 3D effect
						plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, true));

						GenomeAnnotation snpRow = annotationManager.getAnnotation(genome, AnnotationManager.AnnotationType.SNP);
						
						TrackGroup geneGroup = TrackFactory.addGeneTracks(plot,
								createAnnotationDataSource(annotationManager.getAnnotation(
										genome, AnnotationManager.AnnotationType.GENES).getUrl(),
										new GeneParser()),
								createAnnotationDataSource(annotationManager.getAnnotation(
										genome, AnnotationManager.AnnotationType.TRANSCRIPTS).getUrl(),
										new TranscriptParser()),
								createAnnotationDataSource(annotationManager.getAnnotation(
										genome, AnnotationManager.AnnotationType.REFERENCE).getUrl(),
										new SequenceParser()),
								snpRow == null ? null : 
									createAnnotationDataSource(annotationManager.getAnnotation(
											genome, AnnotationManager.AnnotationType.SNP).getUrl(),
											new SNPParser())
								);
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
					switch (track.interpretation.type) {

					case READS:
						if (track.interpretation.summaryDatas.size() == 0) {
							// No precomputed summary data
							TrackFactory.addThickSeparatorTrack(plot);
							treatmentData = createReadDataSource(track.interpretation.primaryData, track.interpretation.indexData, tracks);
							TrackGroup readGroup = TrackFactory.addReadTracks(plot, treatmentData, createReadHandler(file), createAnnotationDataSource(annotationManager.getAnnotation(genome, AnnotationManager.AnnotationType.REFERENCE).getUrl(), new SequenceParser()), track.interpretation.primaryData.getName());
							track.setTrackGroup(readGroup);
						} else { 
							// Has precomputed summary data
							TrackFactory.addThickSeparatorTrack(plot);
							treatmentData = createReadDataSource(track.interpretation.primaryData, track.interpretation.indexData, tracks);
							TrackGroup readGroupWithSummary = TrackFactory.addReadSummaryTracks(plot, treatmentData, createReadHandler(file), createAnnotationDataSource(annotationManager.getAnnotation(genome, AnnotationManager.AnnotationType.REFERENCE).getUrl(), new SequenceParser()), file);
							track.setTrackGroup(readGroupWithSummary);
						}
						break;

					}
				}
			}

			// Add selected peak tracks
			for (Track track : tracks) {
				if (track.checkBox.isSelected()) {
					

					File file = track.interpretation.primaryData == null ? null : Session
							.getSession().getDataManager().getLocalFile(
									track.interpretation.primaryData);
					DataSource peakData;
					switch (track.interpretation.type) {
					case PEAKS:
						TrackFactory.addThickSeparatorTrack(plot);
						TrackFactory.addTitleTrack(plot, file.getName());
						peakData = new ChunkDataSource(file, new BEDParser());
						TrackFactory.addPeakTrack(plot, peakData);
						break;
					case PEAKS_WITH_HEADER:
						TrackFactory.addThickSeparatorTrack(plot);
						TrackFactory.addTitleTrack(plot, file.getName());
						peakData = new ChunkDataSource(file, new HeaderTsvParser());
						TrackFactory.addHeaderPeakTrack(plot, peakData);
						break;
					}
				}
			}

			// End 3D effect
			plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, false));

			// Fill in initial positions if not filled in
			if (locationField.getText().trim().isEmpty()) {
				locationField.setText(DEFAULT_LOCATION);
			}
			if (zoomField.getText().trim().isEmpty()) {
				zoomField.setText(DEFAULT_ZOOM);
			}
			
			// Initialise the plot
			plot.addDataRegionListener(this);

			// Go to correct place (possibly gene name that must be translated)
			updateLocation();
			
			// Remember chromosome
			visibleChromosome = chrBox.getSelectedItem();

			// wrap it in a panel
			chartPanel.setChart(new JFreeChart(plot));
			chartPanel.setCursor(new Cursor(Cursor.HAND_CURSOR));

			// add mouse listeners
			for (View view : plot.getViews()) {
				chartPanel.addMouseListener(view);
				chartPanel.addMouseMotionListener(view);
				chartPanel.addMouseWheelListener(view);
			}

			// put panel on top of card layout
			if (plotPanel.getComponentCount() == 2) {
				plotPanel.remove(1);
			}

			plotPanel.add(chartPanel, PLOTPANEL);
			plotPanel.addComponentListener(this);
			CardLayout cl = (CardLayout) (plotPanel.getLayout());
			cl.show(plotPanel, PLOTPANEL);

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
	 */
	public DataSource createReadDataSource(DataBean data, DataBean indexData, List<Track> tracks)
			throws MicroarrayException, IOException {
		DataSource dataSource = null;

	    // Convert data bean into file
	    File file = data == null ? null : Session.getSession().getDataManager().getLocalFile(data);
	    
	    if (file.getName().contains(".bam-summary")) {
	    	dataSource = new TabixDataSource(file);
	    	
	    } else if (file.getName().contains(".bam") || file.getName().contains(".sam")) {
	    	File indexFile = Session.getSession().getDataManager().getLocalFile(indexData);
	    	dataSource = new SAMDataSource(file, indexFile);
	    	
	    } else {
	    	dataSource = new ChunkDataSource(file, new ElandParser());
	    }
	    
	    return dataSource;
	}

	private boolean isIndexData(DataBean bean) {
		return bean.getName().endsWith(".bai");
	}

    /**
     * Create AreaRequestHandler either for SAM/BAM or ELAND data files.
     * 
     * @param file
     * @return
     */
    public Class<?extends AreaRequestHandler> createReadHandler(File file) {
    	
    	if (file.getName().contains(".bam-summary")) {
    		return TabixHandlerThread.class;
    	}
    	
        if (file.getName().contains(".bam") || file.getName().contains(".sam")) {
            return SAMHandlerThread.class;
        }
        return ChunkTreeHandlerThread.class;
    }

	private ChunkDataSource createAnnotationDataSource(URL url,
			TsvParser fileParser) throws FileNotFoundException, URISyntaxException  {
		if ("file".equals(url.getProtocol())) {
			return new ChunkDataSource(new File(url.toURI()), fileParser);
		} else {
			return new ChunkDataSource(url, fileParser);
		}
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

	public void regionChanged(BpCoordRegion bpRegion) {
		locationField.setText(bpRegion.getMid().toString());
		zoomField.setText("" + bpRegion.getLength());
		this.lastLocation = bpRegion.getMid();
		this.lastZoom = bpRegion.getLength();
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
				interpretations.add(new Interpretation(TrackType.PEAKS, data));

			} else if (data.isContentTypeCompatitible("text/tab")) {
				// peaks (with header in the file)
				interpretations.add(new Interpretation(TrackType.PEAKS_WITH_HEADER, data));

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

			// If gene name was given, search for it and 
			// translate it to coordinates.
			if (!GeneIndexActions.checkIfNumber(locationField.getText())) {

				BpCoordRegion geneLocation = gia.getLocation(locationField.getText().toUpperCase());

				if (geneLocation == null) {
					
					// Move to last known location
					if (lastLocation != null && lastZoom != null) {
						locationField.setText(lastLocation.toString());
						zoomField.setText(lastZoom.toString());
					} else {
						locationField.setText(DEFAULT_LOCATION);
						zoomField.setText(DEFAULT_ZOOM);
					}
					
					// Tell the user 
					application.showDialog("Not found",
							"Gene was not found", null,
							Severity.INFO, true,
							DetailsVisibility.DETAILS_ALWAYS_HIDDEN, null);
					
				} else {
					
					// Update coordinate controls with gene's location
					chrBox.setSelectedItem(new Chromosome(geneLocation.start.chr));
					locationField.setText(Long.toString((geneLocation.end.bp + geneLocation.start.bp) / 2));
					zoomField.setText(Long.toString((geneLocation.end.bp - geneLocation.start.bp) * 2));

				}
			}

			// Check how large the update in location was 
			if (visibleChromosome != null && visibleChromosome != chrBox.getSelectedItem()) {

				// Chromosome changed, redraw everything
				showVisualisation();
				updateVisibilityForTracks();

			} else {

				// Only bp position within chromosome changed, move there
				plot.moveDataBpRegion((Chromosome) chrBox.getSelectedItem(),
						Long.parseLong(locationField.getText()), Long
						.parseLong(zoomField.getText()));

				// Set scale of profile track containing reads information
				this.plot.setReadScale((ReadScale) this.coverageScaleBox.getSelectedItem());
			}		
			
		} catch (Exception e) {
			application.reportException(e);
		}

	}

	public void focusGained(FocusEvent e) {
	}

	public void focusLost(FocusEvent e) {
		// skip
	}

	public void componentHidden(ComponentEvent arg0) {
		// skip
	}

	public void componentMoved(ComponentEvent arg0) {
		// skip
	}

	public void componentResized(ComponentEvent arg0) {
//        showVisualisation();
//        updateVisibilityForTracks();

		this.updateLocation();
		plot.redraw();
	}

	public void componentShown(ComponentEvent arg0) {
		// skip
	}
	
	public ClientApplication getClientApplication() {
		return application;
	}
}
