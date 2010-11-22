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

import org.apache.log4j.Logger;
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
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationContents.Row;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.gbrowser.index.GeneIndexActions;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.IOUtils;

/**
 * Chipster style visualisation for genome browser.
 * 
 * @author Petri Klemel�, Aleksi Kallio
 */
public class GenomeBrowser extends Visualisation implements ActionListener,
		RegionListener, FocusListener, ComponentListener {


	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";

	private static final Logger logger = Logger.getLogger(GenomeBrowser.class);
	
	private static enum TrackType {
		CYTOBANDS(false), 
		GENES(true), 
		TRANSCRIPTS(true), 
		REFERENCE(true),
		PEAKS(true),
		PEAKS_WITH_HEADER(true), 
		READS(true),
		READS_WITH_SUMMARY(true),
		HIDDEN(false);
		
		private boolean isToggleable;

		private TrackType(boolean toggleable) {
			this.isToggleable = toggleable;
		}
	}

	private static class Track {

		TrackType type;
		JCheckBox checkBox;
		String name;
		DataBean userData;
		TrackGroup trackGroup = null;

		public Track(String name, TrackType type) {
			this.name = name;
			this.type = type;
		}

		public Track(String name, TrackType type, DataBean userData) {
			this(name, type);
			this.userData = userData;
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

	private List<DataBean> datas;
	private List<TrackType> interpretations;
	private List<Track> tracks = new LinkedList<Track>();

	private GenomePlot plot;

	private JPanel paramPanel;
	private JPanel settingsPanel = new JPanel();
	private JPanel plotPanel = new JPanel(new CardLayout());

	private JButton gotoButton = new JButton("Go to location");
	private JButton drawButton = new JButton("Draw");

	private JLabel locationLabel = new JLabel("Location (gene or position)");
	private JTextField locationField = new JTextField();

	private JLabel zoomLabel = new JLabel("Zoom");
	private JTextField zoomField = new JTextField(10);
	
	private JLabel chrLabel = new JLabel("Chromosome");
	private JComboBox chrBox = new JComboBox();
	
	private JComboBox genomeBox = new JComboBox();
	
	private JLabel tracksLabel = new JLabel("Datasets");

	private Object lastChromosome;

	private GridBagConstraints settingsGridBagConstraints;
	private AnnotationContents annotationContents;

	private JLabel coverageScaleLabel = new JLabel("Coverage scale");
	private JComboBox coverageScaleBox = new JComboBox();

	private GeneIndexActions gia;

	private boolean initialised;
	
	private Map<JCheckBox, String> trackSwitches = new LinkedHashMap<JCheckBox, String>();
	private Set<JCheckBox> datasetSwitches = new HashSet<JCheckBox>();
	
	
	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);

		// initialize annotations
		this.annotationContents = new AnnotationContents();
		this.annotationContents.initialize();
		
		trackSwitches.put(new JCheckBox("Reads", true), "Reads");
		trackSwitches.put(new JCheckBox("Highlight SNPs", false), "highlightSNP");
		trackSwitches.put(new JCheckBox("Coverage and SNP's", true), "ProfileSNPTrack");
		trackSwitches.put(new JCheckBox("Strand-specific coverage", false), "ProfileTrack");
		trackSwitches.put(new JCheckBox("Quality coverage", false), "QualityCoverageTrack");
		trackSwitches.put(new JCheckBox("Density graph", false), "GelTrack");
//		trackSwitches.put(new JCheckBox("Show reference SNP's", false), "changeSNP"); // TODO re-enable SNP view
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

//		Genome genome = (Genome) genomeBox.getSelectedItem();

		// for now just always add genes and cytobands

//		list available track types for the genome
//		for (Row row : annotationContents.getRows()) {
//			if (genome.equals(row.getGenome())) {
//				TrackType type;
//				if (row.content == AnnotationContents.Content.GENES) {
//					type = TrackType.GENES;
//				} else if (row.content == AnnotationContents.Content.TRANSCRIPTS) {
//					continue; // track not directly supported, skip
//				} else if (row.content == AnnotationContents.Content.CYTOBANDS) {
//					type = TrackType.CYTOBANDS;
//				} else if (row.content == AnnotationContents.Content.REFERENCE) {
//					continue; // track not directly supported, skip
//				} else {
//					continue; // track not supported, skip
//				}
//
//				tracks.add(new Track(row.content.getId(), type));
//			}
//		}

		tracks.add(new Track(AnnotationContents.Content.GENES.getId(), TrackType.GENES));
		tracks.add(new Track(AnnotationContents.Content.CYTOBANDS.getId(), TrackType.CYTOBANDS));
		
		
		for (int i = 0; i < interpretations.size(); i++) {
			TrackType interpretation = interpretations.get(i);
			tracks.add(new Track(datas.get(i).getName(), interpretation, datas
					.get(i)));
		}

		this.settingsGridBagConstraints.gridy++;
		GridBagConstraints c = this.settingsGridBagConstraints;

		
		// draw button
		c.gridy++;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.CENTER;
		drawButton.setEnabled(false);
		settingsPanel.add(drawButton, c);
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1.0;

		
		// add tracks to settings panel
		c.weighty = 0.0;
		
		this.tracksLabel.setEnabled(false);
		settingsPanel.add(tracksLabel, c);
		c.gridy++;
		c.weighty = 1.0;
		c.fill = GridBagConstraints.BOTH;
		c.anchor = GridBagConstraints.NORTHWEST;
		
		JPanel trackPanel = new JPanel();
		trackPanel.setLayout(new BoxLayout(trackPanel, BoxLayout.Y_AXIS));

		for (Track track : tracks) {
			JCheckBox box = new JCheckBox(track.name, true);
			box.setToolTipText(track.name);
			box.setEnabled(false);
			track.checkBox = box;
			if (track.type.isToggleable) {
				trackPanel.add(box);
				datasetSwitches.add(box);
				box.addActionListener(this);
			}
		}
		
		JScrollPane trackPanelScrollPane = new JScrollPane(trackPanel);
		trackPanelScrollPane.setBorder(BorderFactory.createEmptyBorder());
		
		trackPanelScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);		
		settingsPanel.add(trackPanelScrollPane, c);
		c.gridy++;

		
		// options
		c.weighty = 1.0;
		createTracksSwitches();

		// scale options for profile track
		c.fill = GridBagConstraints.HORIZONTAL;
		c.anchor = GridBagConstraints.NORTHWEST;
		c.weighty = 0.0;
		c.gridx = 0;
		c.gridwidth = 5;
		c.gridy++;
		c.insets.set(5, 10, 5, 10);
		coverageScaleLabel.setEnabled(false);
		settingsPanel.add(coverageScaleLabel, c);
		coverageScaleBox = new JComboBox(GenomePlot.ReadScale.values());
		coverageScaleBox.setEnabled(false);
		coverageScaleBox.addActionListener(this);
		c.gridx = 0;
		c.gridwidth = 5;
		c.gridy++;
		settingsPanel.add(coverageScaleBox, c);

	}
	
	public void createTracksSwitches() {
		
		GridBagConstraints c = this.settingsGridBagConstraints;
		
		JPanel menu = new JPanel();
		JScrollPane menuScrollPane = new JScrollPane(menu);
		menuScrollPane.setHorizontalScrollBarPolicy(ScrollPaneConstants.HORIZONTAL_SCROLLBAR_NEVER);
		menuScrollPane.setBorder(BorderFactory.createEmptyBorder());
//		menu.setLayout(new GridLayout(7,1));
		menu.setLayout(new BoxLayout(menu, BoxLayout.Y_AXIS));


		setTrackSwitchesEnabled(false);
		for (JCheckBox trackSwitch : trackSwitches.keySet()) {
			trackSwitch.addActionListener(this);
			menu.add(trackSwitch);
		}
        
		settingsPanel.add(menuScrollPane, c);
		
		c.gridy++;
		c.fill = GridBagConstraints.BOTH;
		c.weighty = 1.0;
		
	}
	
	public JPanel createSettingsPanel() {

		settingsPanel.setLayout(new GridBagLayout());
		settingsPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);

		drawButton.addActionListener(this);

		GridBagConstraints c = new GridBagConstraints();

		c.gridy = 0;
		c.gridx = 0;
		c.insets.set(5, 10, 5, 10);
		c.anchor = GridBagConstraints.NORTHWEST;
		c.fill = GridBagConstraints.HORIZONTAL;
		c.weighty = 0;
		c.weightx = 1.0;
		c.gridx = 0;
		c.gridwidth = 5;

		// genome
		Collection<Genome> genomes = annotationContents.getGenomes();
		c.gridy++;
		settingsPanel.add(new JLabel("Genome"), c);
		
		c.gridy++;
		for (Genome genome : genomes) {
			genomeBox.addItem(genome);
		}
		
		// no selection at startup
		genomeBox.setSelectedItem(null);

		genomeBox.addActionListener(this);
		settingsPanel.add(genomeBox, c);
		c.gridy++;

		// add annotations button
//		JButton button = new JButton("Annotations");
//		button.addActionListener(new ActionListener() {
//
//			@Override
//			public void actionPerformed(ActionEvent arg0) {
//				AnnotationsScreen screen = new AnnotationsScreen(annotationContents);
//				screen.getFrame().pack();
//				screen.getFrame().setVisible(true);
//			}
//		});
//		settingsPanel.add(button, c);
		
		c.gridy++;
		chrLabel.setEnabled(false);
		settingsPanel.add(chrLabel, c);
		c.gridy++;
		chrBox.setEnabled(false);
		settingsPanel.add(chrBox, c);

		// location
		c.gridy++;
		locationLabel.setEnabled(false);
		settingsPanel.add(locationLabel, c);
		c.gridy++;
		locationField.setEnabled(false);
		settingsPanel.add(locationField, c);

		// zoom
		c.gridx = 0;
		c.gridwidth = 5;
		c.gridy++;
		c.insets.set(5, 10, 5, 10);
		zoomLabel.setEnabled(false);
		settingsPanel.add(zoomLabel, c);
		c.gridwidth = 4;
		c.gridy++;
		zoomField.setEnabled(false);
		settingsPanel.add(this.zoomField, c);
		this.zoomField.addFocusListener(this);


		gotoButton.addActionListener(this);
		gotoButton.setEnabled(false);

		this.settingsGridBagConstraints = c;

		return settingsPanel;
	}

	private void fillChromosomeBox() throws IOException {
		TreeSet<String> chromosomes = new TreeSet<String>(); 
		for (int i = 0; i < interpretations.size(); i++) {
			TrackType trackType = interpretations.get(i);
			if (trackType == TrackType.READS || trackType == TrackType.READS_WITH_SUMMARY) {
				DataBean data = datas.get(i);
				InputStream in = null;
				try {
					in  = data.getContentByteStream();
					chromosomes.addAll(SamBamUtils.readChromosomeNames(in));
				} finally {
					IOUtils.closeIfPossible(in);
				}
			}
		}
		for (String chromosome : chromosomes) {
			chrBox.addItem(new Chromosome(chromosome));
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

		if (source == drawButton) {

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
									
									// Make track switches enabled
									setTrackSwitchesEnabled(true);
									
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
			
		} else if (datasetSwitches.contains(source) || source == coverageScaleBox) {
	        showVisualisation();
	        updateVisibilityForTracks();

		} else if (trackSwitches.keySet().contains(source)) {
			updateVisibilityForTracks();
		} 
		
		// genome selected
		else if (source == genomeBox) {
			
			Genome genome = (Genome) genomeBox.getSelectedItem();

			// dialog for downloading annotations if not already local
			if (!annotationContents.hasLocalAnnotations(genome)) {
				annotationContents.openDownloadAnnotationsDialog(genome);
			}

			// enable other settings
			this.drawButton.setEnabled(true);
			this.chrLabel.setEnabled(true);
			this.chrBox.setEnabled(true);
			this.locationLabel.setEnabled(true);
			this.locationField.setEnabled(true);
			this.zoomLabel.setEnabled(true);
			this.zoomField.setEnabled(true);
			
			this.tracksLabel.setEnabled(true);
			for (Track track : tracks) {
				track.checkBox.setEnabled(true);
			}
			
			coverageScaleLabel.setEnabled(true);
			coverageScaleBox.setEnabled(true);
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
		
		this.datas = datas;
		this.interpretations = interpretUserDatas(this.datas);
		
		// List available chromosomes from user data files
		fillChromosomeBox();

		// We can create tracks now that we know the data
		createAvailableTracks(); 

		// Create panel with card layout and put message panel there
		JPanel waitPanel = new JPanel();
		waitPanel.add(new JLabel("Please select settings"));
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
				gia = GeneIndexActions.getInstance(genome, createAnnotationDataSource(annotationContents.getRow(genome, AnnotationContents.Content.GENES).getUrl(),	new GeneParser()));
			} catch (Exception e) {
				logger.warn("could not create gene name index", e);
			}
			
			// Create the plot
			ChartPanel chartPanel = new NonScalableChartPanel();
			this.plot = new GenomePlot(chartPanel, true);
			
			// Set scale of profile track containing reads information
			this.plot.setReadScale((ReadScale) this.coverageScaleBox
					.getSelectedItem());


			// Add selected annotation tracks
			for (Track track : tracks) {
				if (track.checkBox.isSelected()) {
					switch (track.type) {
					case CYTOBANDS:
						TrackFactory.addCytobandTracks(plot,
								createAnnotationDataSource(
										annotationContents.getRow(
												genome, AnnotationContents.Content.CYTOBANDS).getUrl(),
										new CytobandParser()));
						break;
						
					case GENES:
						// Start 3D effect
						plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, true));

						Row snpRow = annotationContents.getRow(genome, AnnotationContents.Content.SNP);
						
						TrackGroup geneGroup = TrackFactory.addGeneTracks(plot,
								createAnnotationDataSource(annotationContents.getRow(
										genome, AnnotationContents.Content.GENES).getUrl(),
										new GeneParser()),
								createAnnotationDataSource(annotationContents.getRow(
										genome, AnnotationContents.Content.TRANSCRIPTS).getUrl(),
										new TranscriptParser()),
								createAnnotationDataSource(annotationContents.getRow(
										genome, AnnotationContents.Content.REFERENCE).getUrl(),
										new SequenceParser()),
								snpRow == null ? null : 
									createAnnotationDataSource(annotationContents.getRow(
											genome, AnnotationContents.Content.SNP).getUrl(),
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


					File file = track.userData == null ? null : Session
							.getSession().getDataManager().getLocalFile(
									track.userData);
					DataSource treatmentData;
					switch (track.type) {

					case READS:
						TrackFactory.addThickSeparatorTrack(plot);
						treatmentData = createReadDataSource(track.userData, tracks);
						TrackGroup readGroup = TrackFactory.addReadTracks(plot, treatmentData, createReadHandler(file), createAnnotationDataSource(annotationContents.getRow(genome, AnnotationContents.Content.REFERENCE).getUrl(), new SequenceParser()), file.getName());
						track.setTrackGroup(readGroup);
						break;

					case READS_WITH_SUMMARY:
						TrackFactory.addThickSeparatorTrack(plot);
						treatmentData = createReadDataSource(track.userData, tracks);
						TrackGroup readGroupWithSummary = TrackFactory.addReadSummaryTracks(plot, treatmentData, createReadHandler(file), createAnnotationDataSource(annotationContents.getRow(genome, AnnotationContents.Content.REFERENCE).getUrl(), new SequenceParser()), file.getName());
						track.setTrackGroup(readGroupWithSummary);
						break;
					}
				}
			}

			// Add selected peak tracks
			for (Track track : tracks) {
				if (track.checkBox.isSelected()) {
					

					File file = track.userData == null ? null : Session
							.getSession().getDataManager().getLocalFile(
									track.userData);
					DataSource peakData;
					switch (track.type) {
					case PEAKS:
						TrackFactory.addThickSeparatorTrack(plot);
						peakData = new ChunkDataSource(file, new BEDParser());
						TrackFactory.addThickSeparatorTrack(plot);
						TrackFactory.addTitleTrack(plot, file.getName());
						TrackFactory.addPeakTrack(plot, peakData);
						break;
					case PEAKS_WITH_HEADER:
						TrackFactory.addThickSeparatorTrack(plot);
						peakData = new ChunkDataSource(file,
								new HeaderTsvParser());
						TrackFactory.addThickSeparatorTrack(plot);
						TrackFactory.addTitleTrack(plot, file.getName());
						TrackFactory.addHeaderPeakTrack(plot, peakData);
						break;
					}
				}
			}

			// End 3D effect
			plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, false));

			// Fill in initial positions if not filled in
			if (locationField.getText().trim().isEmpty()) {
				locationField.setText("1000000");
			}
			if (zoomField.getText().trim().isEmpty()) {
				zoomField.setText("100000");
			}

			// Initialise the plot
			plot.addDataRegionListener(this);

			// remember the chromosome, so we know if it has changed
			lastChromosome = chrBox.getSelectedItem();
			updateLocation();

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
		for (DataBean data : datas) {
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
	public DataSource createReadDataSource(DataBean data, List<Track> tracks)
			throws MicroarrayException, IOException {
		DataSource dataSource = null;

	    // Convert data bean into file
	    File file = data == null ? null : Session.getSession().getDataManager().getLocalFile(data);
	    
	    if (file.getName().contains(".bam-summary")) {
	    	dataSource = new TabixDataSource(file);
	    	
	    } else if (file.getName().contains(".bam") || file.getName().contains(".sam")) {
	    	// Find the index file from the operation
	    	DataBean indexBean = null;
	    	for (Track track : tracks) {
	    		if (track.type == GenomeBrowser.TrackType.HIDDEN) {
	    			DataBean bean = track.userData;
	    			if (isIndexData(bean) && bean.getName().startsWith(data.getName())) {
	    				indexBean = bean;
	    			}
	    		}
	    	}
	    	if (indexBean == null) {
	    		throw new MicroarrayException("Index file not selected for SAM/BAM file " + data.getName());
	    	}
	    	File indexFile = Session.getSession().getDataManager().getLocalFile(indexBean);
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
		// Check if there's non-compatible data  
		for (DataBean data : datas) {
			if (!data.hasTypeTag(MicroarrayModule.TypeTags.ORDERED_GENOMIC_ENTITIES)) {
				return false;
			}
		}
		return true;
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
		gotoButton.setEnabled(false);
	}

	private List<TrackType> interpretUserDatas(List<DataBean> datas) {
		LinkedList<TrackType> interpretations = new LinkedList<TrackType>();

		// try to find interpretation for all selected datas
		for (DataBean data : datas) {

			if (data.isContentTypeCompatitible("text/plain")) {
				interpretations.add(TrackType.READS);

			} else if (data.isContentTypeCompatitible("text/bed")) {
				// peaks
				interpretations.add(TrackType.PEAKS);

			} else if (data.isContentTypeCompatitible("text/tab")) {
				// peaks (with header in the file)
				interpretations.add(TrackType.PEAKS_WITH_HEADER);

			} else if ((data.isContentTypeCompatitible("application/octet-stream")) &&
					(data.getName().contains(".bam-summary"))) {
				interpretations.add(TrackType.READS_WITH_SUMMARY);
				
			} else if ((data.isContentTypeCompatitible("application/octet-stream")) &&
			           (data.getName().endsWith(".bam"))) {
                interpretations.add(TrackType.READS);
                
			} else if ((data.isContentTypeCompatitible("application/octet-stream")) &&
			           (isIndexData(data))) {
				interpretations.add(TrackType.HIDDEN);
                
			} else {
	             throw new RuntimeException("cannot visualise: " + data.getName());
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

		// Chromosome changed - redraw (alternatively we could clean track
		// contents)
		if (lastChromosome != chrBox.getSelectedItem()) {
			showVisualisation();
			updateVisibilityForTracks();
			return;
		}

		// Only position within chromosome changed
		BpCoordRegion geneLocation;
		if (!GeneIndexActions.checkIfNumber(locationField.getText())) {
			if (gia == null) {
				// TODO tell user about it
				return;
			}
			
			
			geneLocation = gia.getLocation(locationField.getText().toUpperCase(), (Chromosome)chrBox.getSelectedItem());

			if (geneLocation == null) {
				application.showDialog("Not found",
						"Gene with such name was not found", null,
						Severity.INFO, false,
						DetailsVisibility.DETAILS_ALWAYS_HIDDEN, null);
			} else {
				chrBox.setSelectedItem(new Chromosome(geneLocation.start.chr));
				plot.moveDataBpRegion((Chromosome) chrBox.getSelectedItem(),
						(geneLocation.end.bp + geneLocation.start.bp) / 2,
						(geneLocation.end.bp - geneLocation.start.bp) * 2);
			}
		} else {
			try {
				plot.moveDataBpRegion((Chromosome) chrBox.getSelectedItem(),
						Long.parseLong(locationField.getText()), Long
								.parseLong(zoomField.getText()));
			} catch (NumberFormatException e) {
				application.reportException(e);
			}
		}

		// Set scale of profile track containing reads information
		this.plot.setReadScale((ReadScale) this.coverageScaleBox
				.getSelectedItem());
	}

	public void focusGained(FocusEvent e) {
		gotoButton.setEnabled(true);
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
