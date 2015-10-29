package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;

import net.miginfocom.swing.MigLayout;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.operation.parameter.SteppedComboBox;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager.AnnotationType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Interpretation.TrackType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.AnnotationTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SampleTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Selectable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;
import fi.csc.microarray.util.LinkUtil;

/**
 * Genome browser settings GUI. The side panel for genome browser to change the settings of the visualisation.
 * 
 * @author klemela
 */
public class GBrowserSettings implements ActionListener, RegionListener {
	
	public enum CoverageType {
		NONE ("none"),
		TOTAL ("total"),
		STRAND ("strand-specific"),
		STRAND_XS ("strand-specific for RNA-seq");
		
		String name;
		
		CoverageType(String name) {
			this.name = name;
		}
		
		public String toString() {
			return name;
		}
		
		public String getId() {
			return "coverage type " + toString();
		}
	}
	
	
	private static final long DEFAULT_VIEWSIZE = 100_000;
	private static final long DEFAULT_LOCATION = 45_000;
	private static final String FULL_WIDTH = "growx";
	private static final String GAPY = "gapy rel";
	
	private Long lastViewsize;
	private boolean initialised;
	
	private JPanel paramPanel;
	private JPanel settingsPanel = new JPanel();

	private JButton goButton = new JButton("Go");

	private JLabel locationLabel = new JLabel("Location (gene or position)");
	private JTextField locationField = new JTextField();

	private static final String VIEW_SIZE_TEXT = "View size: ";
	private JLabel viewsizeLabel = new JLabel(VIEW_SIZE_TEXT);

	private JLabel chrLabel = new JLabel("Chromosome");
	private JComboBox<Chromosome> chrBox;

	private SteppedComboBox genomeBox;

	private JLabel coverageScaleLabel = new JLabel("Coverage scale");
	private JComboBox<ReadScale> coverageScaleBox;

	private JLabel coverageTypeLabel = new JLabel("Coverage type");
	private JComboBox<CoverageType> coverageTypeBox; 

	private Map<String, JCheckBox> trackSwitches = new LinkedHashMap<>();

	private JXHyperlink ensemblLink;
	private JXHyperlink ucscLink;
	private GBrowser browser;
	private Long lastLocation;

	private JTabbedPane tabPane;
	private JScrollPane selectedScrollPane;
	private JScrollPane settingsScrollPane;
	private JScrollPane legendScrollPane;

	public void initialise(GBrowser browser) throws Exception {
		
		this.browser = browser;

		trackSwitches.put("Reads", new JCheckBox("Reads", true));
		trackSwitches.put("highlightSNP", new JCheckBox("Highlight SNPs", true));
		//		trackSwitches.put(new JCheckBox("Quality coverage", false), "QualityCoverageTrack"); // TODO re-enable quality coverage
		trackSwitches.put("DensityGraphTrack", new JCheckBox("Density graph", false));
		trackSwitches.put("RepeatMaskerTrack", new JCheckBox("Low complexity regions", false));
		trackSwitches.put("multimapping", new JCheckBox("Mark multimapping reads", false));
	}
	
	public JPanel getParameterPanel() {

		if (paramPanel == null) {
			paramPanel = new JPanel();
			paramPanel.setLayout(new MigLayout("insets 0, fill"));
			
			settingsPanel.setLayout(new MigLayout("wrap 1, fillx, gapy 0"));
			
			JPanel settings = this.createSettingsPanel();
			JPanel selectionPanel = new SelectionPanel(browser.getSelectionManager()); 
			JPanel legend = new GBrowserLegend(browser);
			
			settingsScrollPane = new JScrollPane(settings);
			selectedScrollPane = new JScrollPane(selectionPanel);
			legendScrollPane = new JScrollPane(legend);

			settingsScrollPane.setBorder(BorderFactory.createEmptyBorder());
			selectedScrollPane.setBorder(BorderFactory.createEmptyBorder());			
			legendScrollPane.setBorder(BorderFactory.createEmptyBorder());
			
			settingsScrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
			

			tabPane = new JTabbedPane();
			tabPane.addTab("Settings", settingsScrollPane);
			tabPane.addTab("Selected", selectedScrollPane);
			tabPane.addTab("Legend", legendScrollPane);
			
			browser.getSelectionManager().addSelectionListener(new BrowserSelectionListener() {				
				@Override
				public void selectionChanged(DataUrl data, Selectable selectable, Object source) {
					//there is nothing to show if the selection was cleared 
					if (selectable != null) {
						tabPane.setSelectedComponent(selectedScrollPane);
					}
				}
			});

			paramPanel.add(tabPane, "top, grow");
		}
		
		return paramPanel;
	}

	public void createOptions() {

		setTrackSwitchesEnabled(false);
		for (JCheckBox trackSwitch : trackSwitches.values()) {
			trackSwitch.addActionListener(this);
			settingsPanel.add(trackSwitch, GAPY);
		}

		// coverage type
		coverageTypeLabel.setEnabled(false);
		settingsPanel.add(coverageTypeLabel, GAPY);

		coverageTypeBox = new JComboBox<CoverageType>(new CoverageType[] {CoverageType.NONE, CoverageType.TOTAL, CoverageType.STRAND, CoverageType.STRAND_XS});
		coverageTypeBox.setSelectedItem(CoverageType.TOTAL);
		coverageTypeBox.setEnabled(false);
		coverageTypeBox.addActionListener(this);
		settingsPanel.add(coverageTypeBox, FULL_WIDTH);

		// coverage scale
		coverageScaleLabel.setEnabled(false);
		settingsPanel.add(coverageScaleLabel, GAPY);

		coverageScaleBox = new JComboBox<ReadScale>(GBrowserPlot.ReadScale.values());
		coverageScaleBox.setSelectedItem(GBrowserPlot.ReadScale.SMALL);
		coverageScaleBox.setEnabled(false);
		coverageScaleBox.addActionListener(this);
		settingsPanel.add(coverageScaleBox, FULL_WIDTH);
	}

	private void createExternalLinks() {
		ensemblLink = LinkUtil.createLink("Ensembl", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {					
				browser.openExternalBrowser(browser.getExternalLinkUrl(AnnotationType.ENSEMBL_BROWSER_URL));
			}
		});


		ucscLink = LinkUtil.createLink("UCSC", new AbstractAction() {
			@Override
			public void actionPerformed(ActionEvent arg0) {			
				browser.openExternalBrowser(browser.getExternalLinkUrl(AnnotationType.UCSC_BROWSER_URL));
			}
		});

		ensemblLink.setEnabled(false);
		ucscLink.setEnabled(false);
		
		/* there shouldn't be need for two panels, because layout constraint
		 * "nogrid" separates columns of the two rows, but it seems to 
		 * disable "gap 0" settings
		 */
		
		JPanel row1 = new JPanel(new MigLayout("insets 0, gap 0"));
		JPanel row2 = new JPanel(new MigLayout("insets 0, gap 0"));
					
		row1.add(new JLabel("View this region in "));
		row1.add(ensemblLink);
		row1.add(new JLabel(" or "), "wrap");
		row2.add(ucscLink);
		row2.add(new JLabel(" genome browser."));
		
		settingsPanel.add(row1, GAPY);
		settingsPanel.add(row2, GAPY);
	}


	private void setExternalLinksEnabled() {

		ensemblLink.setEnabled(browser.getExternalLinkUrl(AnnotationType.ENSEMBL_BROWSER_URL).length() > 0);
		ucscLink.setEnabled(browser.getExternalLinkUrl(AnnotationType.UCSC_BROWSER_URL).length() > 0);
	}

	public JPanel createSettingsPanel() {

		settingsPanel.add(createTitle("Genome"), "gaptop unrelated");
		createGenomeSettings();
		settingsPanel.add(createTitle("Location"), "gaptop unrelated");
		createLocationSettings();
		settingsPanel.add(createTitle("Options"), "gaptop unrelated");
		createOptions();
		settingsPanel.add(createTitle("External links"), "gaptop unrelated");
		createExternalLinks();

		return settingsPanel;
	}


	private void createGenomeSettings() {

		Object[] genomes = browser.getAnnotationManager().getGenomes().toArray();
		genomeBox = new SteppedComboBox(genomes);

		genomeBox.setPopupWidth(500);
		genomeBox.setBackground(Color.white);

		// no selection at startup
		genomeBox.setSelectedItem(null);
		genomeBox.addActionListener(this);

		// genomeBox works properly only when min:pref:max sizes are set 
		settingsPanel.add(genomeBox, "w 100:200:300, growx, " + GAPY);
	}

	private void createLocationSettings() {

		// chromosome
		chrLabel.setEnabled(false);
		settingsPanel.add(chrLabel, GAPY);
		chrBox = new JComboBox<Chromosome>();
		chrBox.setEnabled(false);
		settingsPanel.add(chrBox, FULL_WIDTH);

		// location
		locationLabel.setEnabled(false);
		settingsPanel.add(locationLabel, GAPY);
		locationField.setEnabled(false);
		locationField.addActionListener(this);
		settingsPanel.add(locationField, FULL_WIDTH);

		// viewsize
		viewsizeLabel.setEnabled(false);
		setCoordinateFields(null, DEFAULT_VIEWSIZE);
		lastViewsize = DEFAULT_VIEWSIZE;
		settingsPanel.add(viewsizeLabel, GAPY);		

		// go button
		goButton.setEnabled(false);
		goButton.addActionListener(this);
		settingsPanel.add(goButton, FULL_WIDTH + ", " + GAPY);
	}

	protected void fillChromosomeBox() throws IOException, UnsortedDataException, URISyntaxException, GBrowserException {

		LinkedList<Chromosome> chromosomes = browser.getChromosomeNames();
		
		// Fill in the box
		for (Chromosome chromosome : chromosomes) {
			chrBox.addItem(chromosome);
		}
	}
	
	volatile boolean failed = false;
	
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

				browser.runBlockingTask("initialising genome browser", new Runnable() {
					@Override
					public void run() {
						try {					

							// Update UI in Event Dispatch Thread
							SwingUtilities.invokeAndWait(new Runnable() {
								@Override
								public void run() {
									try {
										// Show visualisation
										browser.showVisualisation();
									} catch (Exception e) {
										browser.reportException(e);
										failed = true;
									}
								}
							});
							
							// Update UI in Event Dispatch Thread, update location only after this block task 
							// has quit to be able to show gene search blocking task. It isn't critical if the timing
							// fails, search will still work, but the fancy white glass pane won't show.
							if (!failed) {
								SwingUtilities.invokeLater(new Runnable() {
									@Override
									public void run() {

										try {
											processLocationPanelInput();
											setExternalLinksEnabled();

											initialised = true;

										} catch (Exception e) {
											browser.reportException(e);
										}
									}
								});
							}

						} catch (Exception e) {
							throw new RuntimeException(e);
						}
					}
				});

			} else {

				// Move to correct location
				processLocationPanelInput();

				this.setExternalLinksEnabled();
			}

		} else if (source == coverageScaleBox && this.initialised) {

			browser.updateCoverageScale();

		} else if ((trackSwitches.values().contains(source) || source == coverageTypeBox) && this.initialised) {
			updateVisibilityForTracks();
		} 

		// genome selected
		else if (source == genomeBox) {

			Genome genome = (Genome) genomeBox.getSelectedItem();

			// dialog for downloading annotations if not already local
			if (!browser.getAnnotationManager().hasLocalAnnotations(genome)) {
				browser.openDownloadAnnotationsDialog(genome);
			}

			// enable other settings
			this.goButton.setEnabled(true);
			this.chrLabel.setEnabled(true);
			this.chrBox.setEnabled(true);
			this.locationLabel.setEnabled(true);
			this.locationField.setEnabled(true);
			this.viewsizeLabel.setEnabled(true);

			coverageTypeLabel.setEnabled(true);
			coverageTypeBox.setEnabled(true);

			coverageScaleLabel.setEnabled(true);
			coverageScaleBox.setEnabled(true);

			this.setTrackSwitchesEnabled(true);
		}
	}
	
	private void setTrackSwitchesEnabled(boolean enabled) {
		for (JCheckBox trackSwitch : trackSwitches.values()) {
			trackSwitch.setEnabled(enabled);
			
			if (enabled && "RepeatMaskerTrack".equals(trackSwitches.get(trackSwitch))) {
				trackSwitch.setEnabled(browser.getAnnotationUrl(getGenome(), AnnotationType.REPEAT) != null);
			}
		}
	}

	/**
	 * Null keeps the existing content.
	 * 
	 * @param location
	 * @param viewsize
	 */
	public void setCoordinateFields(Long location, Long viewsize) {
		if (location != null) {
			locationField.setText(location.toString());
		}

		if (viewsize != null) {
			if (viewsize > 1000000) {
				viewsizeLabel.setText(VIEW_SIZE_TEXT + Math.round(((float)viewsize) / 1000000f) + " Mb");
			} else if (viewsize > 1000) {
				viewsizeLabel.setText(VIEW_SIZE_TEXT + Math.round(((float)viewsize) / 1000f) + " kb");
			} else {
				viewsizeLabel.setText(VIEW_SIZE_TEXT + viewsize + "");
			}
		}
	}
	
	/**
	 * Changes the visibility of some tracks. This is useful when track is hidden to free some 
	 * screen estate or drawing performance, but the data is still kept in memory (because it's 
	 * needed elsewhere or we just don't care about it).
	 * 
	 * Currently this is used to change between different views of user's bam-files.
	 * 
	 * For more radical changes use updateTracks(), which removes the data layers as well.
	 * 
	 */
	public void updateVisibilityForTracks() {
		for (TrackGroup trackGroup : browser.getPlot().getDataView().getTrackGroups()) {
			
			if (trackGroup instanceof AnnotationTrackGroup) {
				AnnotationTrackGroup annotations = (AnnotationTrackGroup) trackGroup;
				annotations.setRepeatVisible(trackSwitches.get("RepeatMaskerTrack").isSelected());
			}

			if (trackGroup instanceof SampleTrackGroup) {
				SampleTrackGroup samples = (SampleTrackGroup)trackGroup;
				samples.setCoverageType((CoverageType) coverageTypeBox.getSelectedItem());
				samples.setReadsVisible(trackSwitches.get("Reads").isSelected());
				samples.setHighlightSnp(trackSwitches.get("highlightSNP").isSelected());
				samples.setMarkMultimappingReads(trackSwitches.get("multimapping").isSelected());
				samples.setDensityGraphVisible(trackSwitches.get("DensityGraphTrack").isSelected());
			}		
		}
		this.browser.getPlot().getDataView().reloadData();
	}

	public String getGoButtonText() {
		return goButton.getText();
	}

	public ReadScale getCoverageScale() {
		return (ReadScale)coverageScaleBox.getSelectedItem();
	}

	public Genome getGenome() {
		return (Genome) genomeBox.getSelectedItem();
	}

	public Chromosome getChromosome() {
		return (Chromosome)chrBox.getSelectedItem();
	}

	public void processLocationPanelInput() {
		
		// Fill in initial position if not filled in
		if (locationField.getText().trim().isEmpty()) {

			setCoordinateFields(DEFAULT_LOCATION, null);
		}
		
		if (!locationField.getText().isEmpty()) {
			if (!GeneIndexActions.checkIfNumber(locationField.getText())) {

				// If gene name was given, search for it
				browser.requestGeneSearch(locationField.getText());
				setCoordinateFields(getLocation(), getViewSize());
			} else {
				browser.setLocation(getChromosome(), getLocationStart(), getLocationEnd());
			}
		}	
	}


	public Long getLocation() {
		try {
			return Long.parseLong(locationField.getText());
		} catch (NumberFormatException e) {
			if (lastLocation != null) {
				return lastLocation;
			} else {
				return DEFAULT_LOCATION;
			}
		}
	}
	
	private Long getLocationStart() {
		return getLocation() - getViewSize() / 2;
	}
	
	private Long getLocationEnd() {
		return getLocation() + getViewSize() / 2;
	}

	public Long getViewSize() {
		if (lastViewsize != null) {
			return lastViewsize;
		} else {
			return DEFAULT_VIEWSIZE;
		}
	}

	public boolean setChromosome(Chromosome chr) {
		
		chrBox.setSelectedItem(chr);

		return chrBox.getSelectedItem().equals(chr);
	}
	

	public void regionChanged(Region bpRegion) {
		setCoordinateFields(bpRegion.getMid(), bpRegion.getLength());
		this.lastLocation = bpRegion.getMid();
		this.lastViewsize = bpRegion.getLength();
	}
	
	public static JLabel createTitle(String title) {
		JLabel label = new JLabel(title);
		label.setFont(UIManager.getFont("TitledBorder.font"));
		label.setForeground(UIManager.getColor("TitledBorder.titleColor"));
		return label;
	}

	@SuppressWarnings("unchecked")
	public void updateInterpretations() throws IOException, UnsortedDataException, URISyntaxException, GBrowserException {
		
		for (Interpretation interpretation : browser.getInterpretations()) {
			if (interpretation.getType() == TrackType.REFERENCE) {
				Genome ownGenome = new Genome(interpretation.getName(), "");
				genomeBox.addItem(ownGenome);
			}
		}
		
		fillChromosomeBox();		
	}

	public void updateTracks() {
		updateVisibilityForTracks();
	}
}

