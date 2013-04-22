package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;

import javax.swing.AbstractAction;
import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.ScrollPaneConstants;
import javax.swing.SwingUtilities;
import javax.swing.border.TitledBorder;

import org.jdesktop.swingx.JXHyperlink;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.TrackDefinition;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.AnnotationType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadTrackGroup;
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
		STRAND ("strand-specific");
		
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
	
	
	private static final long DEFAULT_VIEWSIZE = 100000;
	private static final long DEFAULT_LOCATION = 1000000;
	
	private Long lastViewsize;
	private boolean initialised;
	
	private JPanel paramPanel;
	private JPanel settingsPanel = new JPanel();
	private JPanel genomePanel;
	private JPanel locationPanel;
	private JPanel datasetsPanel;
	private JPanel datasetSwitchesPanel;
	private JPanel optionsPanel;
	private JPanel linksPanel;

	private JButton goButton = new JButton("Go");

	private JLabel locationLabel = new JLabel("Location (gene or position)");
	private JTextField locationField = new JTextField();

	private JLabel viewsizeLabel = new JLabel("View size");
	private JTextField viewsizeField = new JTextField();

	private JLabel chrLabel = new JLabel("Chromosome");
	private JComboBox<Chromosome> chrBox;

	private JComboBox<Genome> genomeBox;

	private JLabel coverageScaleLabel = new JLabel("Coverage scale");
	private JComboBox<ReadScale> coverageScaleBox;

	private JLabel coverageTypeLabel = new JLabel("Coverage type");
	private JComboBox<CoverageType> coverageTypeBox; 

	private Map<JCheckBox, String> trackSwitches = new LinkedHashMap<JCheckBox, String>();
	private Set<JCheckBox> datasetSwitches = new HashSet<JCheckBox>();

	private JCheckBox showFullHeightBox;

	private JXHyperlink ensemblLink;
	private JXHyperlink ucscLink;
	private GBrowser browser;
	private Long lastLocation;
	private JCheckBox cacheBox;

	public void initialise(GBrowser browser) throws Exception {
		
		this.browser = browser;

		trackSwitches.put(new JCheckBox("Reads", true), "Reads");
		trackSwitches.put(new JCheckBox("Highlight SNPs", true), "highlightSNP");
		//		trackSwitches.put(new JCheckBox("Coverage and SNPs", true), "ProfileSNPTrack");
		//		trackSwitches.put(new JCheckBox("Strand-specific coverage", false), "ProfileTrack");

		//		trackSwitches.put(new JCheckBox("Quality coverage", false), "QualityCoverageTrack"); // TODO re-enable quality coverage
		trackSwitches.put(new JCheckBox("Density graph", false), "GelTrack");
		trackSwitches.put(new JCheckBox("Low complexity regions", false), "RepeatMaskerTrack");
		//		trackSwitches.put(new JCheckBox("Known SNP's", false), "changeSNP"); // TODO re-enable dbSNP view
	}
	
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
			tabPane.addTab("Legend", new GBrowserLegend(browser));

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

	private JPanel getDatasetsPanel() {
		if (datasetsPanel == null) { 
			datasetsPanel = new JPanel(new GridBagLayout());
			datasetsPanel.setBorder(createSettingsPanelSubPanelBorder("Datasets"));

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

	public JPanel getOptionsPanel() {
		if (this.optionsPanel == null) {
			optionsPanel = new JPanel(new GridBagLayout());
			optionsPanel.setBorder(createSettingsPanelSubPanelBorder("Options"));

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
			coverageTypeBox = new JComboBox<CoverageType>(new CoverageType[] {CoverageType.NONE, CoverageType.TOTAL, CoverageType.STRAND});
			coverageTypeBox.setSelectedItem(CoverageType.TOTAL);
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
			coverageScaleBox = new JComboBox<ReadScale>(GBrowserPlot.ReadScale.values());
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
			linksPanel.setBorder(createSettingsPanelSubPanelBorder("External links"));

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

		ensemblLink.setEnabled(browser.getExternalLinkUrl(AnnotationType.ENSEMBL_BROWSER_URL).length() > 0);
		ucscLink.setEnabled(browser.getExternalLinkUrl(AnnotationType.UCSC_BROWSER_URL).length() > 0);
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
			genomePanel.setBorder(createSettingsPanelSubPanelBorder("Genome"));

			GridBagConstraints c = new GridBagConstraints();
			c.gridy = 0;
			c.gridx = 0;
			c.insets.set(5, 0, 5, 0);
			c.anchor = GridBagConstraints.NORTHWEST;
			c.fill = GridBagConstraints.HORIZONTAL;
			c.weighty = 0;
			c.weightx = 1.0;
			c.gridx = 0;

			genomeBox = new JComboBox<Genome>();
			
			// genome
			Collection<Genome> genomes = browser.getAnnotationManager().getGenomes();
			for (Genome genome : genomes) {
				genomeBox.addItem(genome);
			}

			// no selection at startup
			genomeBox.setSelectedItem(null);
			genomeBox.addActionListener(this);

			genomePanel.add(genomeBox, c);
			
			c.gridy++;
			
			cacheBox = new JCheckBox("Download workflow data");
			cacheBox.setToolTipText("Download workflow data to create a temporary local copy of all user data files. " +
					"The genome browser works faster with the local data, but it will take some time to download the " +
					"files when the visualisation is started.");
			genomePanel.add(cacheBox, c);
		}

		return genomePanel;
	}

	private JPanel getLocationPanel() {
		if (this.locationPanel == null) {

			locationPanel = new JPanel(new GridBagLayout());
			locationPanel.setBorder(createSettingsPanelSubPanelBorder("Location"));

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
			chrBox = new JComboBox<Chromosome>();
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

	protected void fillChromosomeBox() throws IOException {

		LinkedList<Chromosome> chromosomes = browser.getChromosomeNames();
		
		// Fill in the box
		for (Chromosome chromosome : chromosomes) {
			chrBox.addItem(chromosome);
		}
	}

	public void updateDatasetSwitches() {
		datasetSwitches.clear();

		datasetSwitchesPanel.removeAll();

		for (TrackDefinition track : browser.getTracks()) {
			JCheckBox box = new JCheckBox(track.name, true);
			box.setToolTipText(track.name);
			box.setEnabled(false);
			track.checkBox = box;
			if (track.interpretation.getType().isToggleable) {
				datasetSwitchesPanel.add(box);
				datasetSwitches.add(box);
				box.addActionListener(this);
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
			this.cacheBox.setEnabled(false);
			if (!initialised) {

				browser.runBlockingTask("initialising genome browser", new Runnable() {
					@Override
					public void run() {
						try {
							
							if (cacheBox.isSelected()) {
								// Create a local random access copy of all files in background thread
								browser.initialiseUserDatas();
							}

							// Update UI in Event Dispatch Thread
							SwingUtilities.invokeAndWait(new Runnable() {
								@Override
								public void run() {
									try {
										// Show visualisation
										browser.showVisualisation();
									} catch (Exception e) {
										browser.reportException(e);
									}
								}
							});
							
							// Update UI in Event Dispatch Thread, update location only after this block task 
							// has quit to be able to show gene search blocking task. It isn't critical if the timing
							// fails, search will still work, but the fancy white glass pane won't show.
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

		} else if (datasetSwitches.contains(source) && this.initialised) {

			showFullHeightBox.setSelected(false);
			browser.setFullHeight(false);

			browser.updateData();

		} else if (source == coverageScaleBox && this.initialised) {

			browser.updateCoverageScale();

		} else if ((trackSwitches.keySet().contains(source) || source == coverageTypeBox) && this.initialised) {
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
			this.viewsizeField.setEnabled(true);
			this.showFullHeightBox.setEnabled(true);

			for (TrackDefinition track : browser.getTracks()) {
				track.checkBox.setEnabled(true);
			}

			coverageTypeLabel.setEnabled(true);
			coverageTypeBox.setEnabled(true);

			coverageScaleLabel.setEnabled(true);
			coverageScaleBox.setEnabled(true);

			this.setTrackSwitchesEnabled(true);

		} else if (source == showFullHeightBox && this.initialised) {

			browser.setFullHeight(showFullHeightBox.isSelected());
		}
	}
	
	private void setTrackSwitchesEnabled(boolean enabled) {
		for (JCheckBox trackSwitch : trackSwitches.keySet()) {
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
				viewsizeField.setText(Math.round(((float)viewsize) / 1000000f) + " Mb");
			} else if (viewsize > 1000) {
				viewsizeField.setText(Math.round(((float)viewsize) / 1000f) + " kb");
			} else {
				viewsizeField.setText(viewsize + "");
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
		for (TrackDefinition track : browser.getTracks()) {
			if (track.trackGroup != null) {
				for (JCheckBox trackSwitch : trackSwitches.keySet()) {
					track.trackGroup.showOrHide(trackSwitches.get(trackSwitch), trackSwitch.isSelected());
				}
			
				if (track.trackGroup instanceof ReadTrackGroup) {										
					if (coverageTypeBox.getSelectedItem() == CoverageType.NONE) {
						track.trackGroup.showOrHide("Coverage", false);
						track.trackGroup.showOrHide("CoverageEstimate", false);
					} else 	if (coverageTypeBox.getSelectedItem() == CoverageType.TOTAL) {
						track.trackGroup.showOrHide("Coverage", true);
						track.trackGroup.showOrHide("CoverageEstimate", true);	
						track.trackGroup.showOrHide("StrandSpecificCoverageType", false);
					} else 	if (coverageTypeBox.getSelectedItem() == CoverageType.STRAND) {
						track.trackGroup.showOrHide("Coverage", true);
						track.trackGroup.showOrHide("CoverageEstimate", true);	
						track.trackGroup.showOrHide("StrandSpecificCoverageType", true);
					}
				}
			}
		}
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

	public boolean isFullHeight() {
		return showFullHeightBox.isSelected();
	}

	public void processLocationPanelInput() {
		
		// Fill in initial position if not filled in
		if (locationField.getText().trim().isEmpty()) {

			setCoordinateFields(DEFAULT_LOCATION, null);
		}

		if (viewsizeField.getText().trim().isEmpty()) {

			setCoordinateFields(null, DEFAULT_VIEWSIZE);
			lastViewsize = DEFAULT_VIEWSIZE;
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
	
	public static TitledBorder createSettingsPanelSubPanelBorder(String title) {
		return BorderFactory.createTitledBorder(BorderFactory.createMatteBorder(1, 0, 0, 0, Color.lightGray), title);
	}

	public void updateInterpretations() throws IOException {
		fillChromosomeBox();		
	}

	public void updateTracks() {
		updateVisibilityForTracks();
	}
}
