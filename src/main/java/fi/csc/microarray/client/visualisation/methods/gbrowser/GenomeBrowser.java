package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.CardLayout;
import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.visualisation.NonScalableChartPanel;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomePlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.BedTabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GeneSearchHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GtfTabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixSummaryHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParserWithCoordinateConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.HeaderTsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.VcfParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.index.GeneIndexActions;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.AnnotationType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.GenomeAnnotation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.ReadTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;

public class GenomeBrowser implements RegionListener, ComponentListener {
	
	public static enum TrackType {
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
	
	public static class Interpretation {

		public TrackType type;
		public List<URL> summaryDatas = new LinkedList<URL>();
		public URL primaryData;
		public String primaryDataName;
		public URL indexData;

		public Interpretation(TrackType type, URL primaryData, String primaryDataName) {
			this.type = type;
			this.primaryData = primaryData;
			this.primaryDataName = primaryDataName;
		}
				
		protected InputStream getPrimaryDataInputStream() {
			
			/*
			 * 						DataBean data = interpretation.primaryData;
				InputStream in = null;
				try {
					in  = data.getContentByteStream();
					chromosomeNames.addAll(SamBamUtils.readChromosomeNames(in));
				} finally {
					IOUtils.closeIfPossible(in);
				}
			 */
			return null;
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
	
	private static final long DEFAULT_VIEWSIZE = 100000;
	private static final long DEFAULT_LOCATION = 1000000;
	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";
	
	private List<Track> tracks = new LinkedList<Track>();

	private GenomePlot plot;

	private JPanel plotPanel = new JPanel(new CardLayout());

	private AnnotationManager annotationManager;

	private GeneIndexActions gia;

	private boolean initialised;

	private Long lastViewsize;
	private JScrollPane verticalScroller;
	private ViewLimiter viewLimiter;
	protected boolean geneSearchDone;
	
	private GenomeBrowserSettings settings;
	
	private List<Interpretation> interpretations;

	public void initialise() throws Exception {

		// initialize annotations
		this.annotationManager = new AnnotationManager();
		this.annotationManager.initialize();

		settings.initialise();
	}



	private void createAvailableTracks() {

		// for now just always add genes and cytobands
		tracks.add(new Track(AnnotationManager.AnnotationType.GTF_TABIX.getId(), new Interpretation(TrackType.GENES, null, null)));
		tracks.add(new Track(AnnotationManager.AnnotationType.CYTOBANDS.getId(), new Interpretation(TrackType.CYTOBANDS, null, null)));


		for (int i = 0; i < interpretations.size(); i++) {
			Interpretation interpretation = interpretations.get(i);
			tracks.add(new Track(interpretation.primaryDataName, interpretation));
		}

		// update the dataset switches in the settings panel
		settings.updateDatasetSwitches();
	}

	private void setFullHeight(boolean fullHeight) {

		if (fullHeight) {
			verticalScroller.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		} else {
			verticalScroller.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_NEVER);
		}

		plot.setFullHeight(fullHeight);		
	}

	public JComponent getVisualisation(List<Interpretation> interpretations) {
		
		this.interpretations = interpretations;  

		// We can create tracks now that we know the data
		this.tracks.clear();
		createAvailableTracks(); 

		// Create panel with card layout and put message panel there
		JPanel waitPanel = new JPanel(new GridBagLayout());
		GridBagConstraints c = new GridBagConstraints();

		waitPanel.add(new JLabel("<html><p style=\"font-size: larger\">Please select genome and click " + settings.getGoButtonText() + "</p></html>"), c);
		plotPanel.add(waitPanel, WAITPANEL);

		return plotPanel;
	}
	
	private void updateCoverageScale() {
		// Set scale of profile track containing reads information
		this.plot.setReadScale(settings.getCoverageScale());
	}

	private Genome getGenome() {
		return settings.getGenome();
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

					URL cytobandUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.CYTOBANDS);

					try {
						
						if (cytobandUrl != null) {
							CytobandDataSource cytobandDataSource;
							cytobandDataSource = new CytobandDataSource(cytobandUrl);

							TrackFactory.addCytobandTracks(plot, cytobandDataSource);

							this.viewLimiter = new ViewLimiter(plot.getOverviewView().getQueueManager(), 
									cytobandDataSource, plot.getOverviewView());
							this.plot.getDataView().setViewLimiter(viewLimiter);
							this.plot.getOverviewView().setViewLimiter(viewLimiter);
						}

					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					}

					break;

				case GENES:
					// Start 3D effect
					plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, true));

					URL gtfUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.GTF_TABIX);

					URL gtfIndexUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.GTF_TABIX_INDEX);

					URL repeatUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REPEAT);

					URL repeatIndexUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REPEAT_INDEX);

					TabixDataSource gtfDataSource = null;
					TabixDataSource repeatDataSource = null;

					try {
						if (gtfUrl != null && gtfIndexUrl != null) {
							gtfDataSource = new TabixDataSource(gtfUrl, gtfIndexUrl, GtfTabixHandlerThread.class);
						}

						if (repeatUrl != null && repeatIndexUrl != null) {
							repeatDataSource = new TabixDataSource(repeatUrl, repeatIndexUrl, BedTabixHandlerThread.class);
						}

						//Show ruler track even if there are now data sources
						TrackGroup geneGroup = TrackFactory.addGeneTracks(plot, gtfDataSource, repeatDataSource);
						track.setTrackGroup(geneGroup);

					} catch (URISyntaxException e) {
						reportException(e);
					} catch (IOException e) {
						reportException(e);
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

						URL fastaUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REFERENCE);
						URL fastaIndexUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REFERENCE_INDEX);

						IndexedFastaDataSource refSeqDataSource = null;
						
						if (fastaUrl != null && fastaIndexUrl != null) {
							refSeqDataSource = new IndexedFastaDataSource(fastaUrl, fastaIndexUrl);
						}

						if (track.interpretation.summaryDatas.size() == 0) {
							// No precomputed summary data
							TrackFactory.addThickSeparatorTrack(plot);
							treatmentData = createReadDataSource(track.interpretation.primaryData, track.interpretation.indexData, tracks);

							TrackGroup readGroup = TrackFactory.addReadTracks(
									plot, treatmentData, 
									refSeqDataSource, 
									track.interpretation.primaryDataName);

							track.setTrackGroup(readGroup);

						} else { 
							// Has precomputed summary data
							TrackFactory.addThickSeparatorTrack(plot);
							treatmentData = createReadDataSource(track.interpretation.primaryData, track.interpretation.indexData, tracks);
							TrackGroup readGroupWithSummary = TrackFactory.addReadSummaryTracks(
									plot, treatmentData, refSeqDataSource, 
									track.interpretation.primaryDataName, new TabixDataSource(file.toURI().toURL(), null, TabixSummaryHandlerThread.class));
							track.setTrackGroup(readGroupWithSummary);
						}
					}
				} catch (IOException e) {
					reportException(e);
				} catch (URISyntaxException e) {
					reportException(e);
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
						reportException(e);
					}
				}

				DataSource regionData;
				switch (track.interpretation.type) {
				case REGIONS:
					TrackFactory.addThickSeparatorTrack(plot);
					TrackFactory.addTitleTrack(plot, track.interpretation.primaryDataName);

					try {
						regionData = new ChunkDataSource(fileUrl, new BEDParserWithCoordinateConversion(), ChunkTreeHandlerThread.class);
						((ChunkDataSource)regionData).checkSorting();
						TrackFactory.addPeakTrack(plot, regionData);

					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					} catch (IOException e) {
						reportException(e);
					} catch (UnsortedDataException e) {
						showDialog("Unsorted data", e.getMessage(), null, true, true);
					}
					break;
				case REGIONS_WITH_HEADER:
					TrackFactory.addThickSeparatorTrack(plot);
					TrackFactory.addTitleTrack(plot, track.interpretation.primaryDataName);

					try {
						regionData = new ChunkDataSource(fileUrl, new HeaderTsvParser(), ChunkTreeHandlerThread.class);
						TrackFactory.addPeakTrack(plot, regionData);

					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					}
					break;
				case VCF:
					TrackFactory.addThickSeparatorTrack(plot);
					TrackFactory.addTitleTrack(plot, track.interpretation.primaryDataName);

					try {
						regionData = new ChunkDataSource(fileUrl, new VcfParser(), ChunkTreeHandlerThread.class);
						TrackFactory.addPeakTrack(plot, regionData);

					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					}
					break;

				}
			}
		}

		// End 3D effect
		plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, false));
	}

	private URL getAnnotationUrl(Genome genome, AnnotationManager.AnnotationType type) {
		GenomeAnnotation annotation = annotationManager.getAnnotation(
				genome, type);
		if (annotation != null) {
			return annotation.getUrl();					
		} else {
			return null;
		}
	}

	protected void showVisualisation() {

		//Clean old data layers
		if (plot != null) {
			plot.clean();
		}

		// Create the chart panel with tooltip support				
		TooltipAugmentedChartPanel chartPanel = new TooltipAugmentedChartPanel();
		this.plot = new GenomePlot(chartPanel, true);
		
		//FIXME remove Chipster dependency
		((NonScalableChartPanel)chartPanel).setGenomePlot(plot);

		//Set location to plot to avoid trouble in track initialization. 
		//Can't do this with updateLocation, because it would lose gene search when the 
		//tracks clear all data layers
		plot.getDataView().setBpRegion(new RegionDouble(
				DEFAULT_LOCATION - DEFAULT_VIEWSIZE / 2.0, DEFAULT_LOCATION + DEFAULT_VIEWSIZE / 2.0, 
				settings.getChromosome()), true);
		
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

		// Put panel on top of card layout
		if (plotPanel.getComponentCount() == 2) {
			plotPanel.remove(1);
		}

		verticalScroller = new JScrollPane(chartPanel);
		verticalScroller.getVerticalScrollBar().setUnitIncrement(30);

		setFullHeight(settings.isFullHeight());

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
				reportException(e);
			}
		}
		return gia;
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
		this.plot.setReadScale(settings.getCoverageScale());
	}



	private void requestGeneSearch() {

		runBlockingTask("searching gene", new Runnable() {

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

							showDialog("Search failed",
									"Unexpected error happened in the search. Please inform the developers if the problem persists.", null,
									true,
									false, null);
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
						showDialog("Different chromosome", 
								"Searched gene was found from chromosome " + resultChr + " but there is no data for that chromosome", "" + geneLocation, 
								true, false);
					}
				}
			}
		});
	}
	
	public void setLocation(Chromosome chr, Long start, Long end) {

		// Move to selected region 
		chrBox.setSelectedItem(chr));
		if (end == null) {
			end = start;
		}
		setCoordinateFields((end + start) / 2, (end - start) * 2);

		// Update
		updateLocation();
	}
	
	public void removeVisualisation() {

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
	}
	

	public void regionChanged(Region bpRegion) {
		settings.setCoordinateFields(bpRegion.getMid(), bpRegion.getLength());
		this.lastViewsize = bpRegion.getLength();
	}
	
	protected void getChromosomeNames() throws IOException {

		// Gather all chromosome names from all indexed datasets (SAM/BAM)
		TreeSet<String> chromosomeNames = new TreeSet<String>(); 
		for (Interpretation interpretation : interpretations) {
			if (interpretation.type == TrackType.READS) {
				InputStream in;
				try {
					in  = interpretation.getPrimaryDataInputStream();
					chromosomeNames.addAll(SamBamUtils.readChromosomeNames(in));
				} finally {
					IOUtils.closeIfPossible(in);
				}
			}
		}

		// If we still don't have names, go through non-indexed datasets
		if (chromosomeNames.isEmpty()) {
			for (Interpretation interpretation : browser.getInterpretations()) {
				if (interpretation.type == TrackType.REGIONS) {
					DataBean data = interpretation.primaryData;
					File file = Session.getSession().getDataManager().getLocalFile(data);
					List<RegionContent> rows = null;
					try {
						//FIXME remove Chipster dependency 
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
	
	@Override
	public void componentShown(ComponentEvent arg0) {
		// Ignore
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

	public List<Interpretation> getInterpretations() {
		return interpretations;
	}
	
	protected String getExternalLinkUrl(AnnotationType browser) {
		settings.getGenome();
		URL url = annotationManager.getAnnotation(settings.getGenome(), browser).getUrl();

		if (url != null) {
			String stringUrl = url.toString();
			Region region = plot.getDataView().getBpRegion();
			stringUrl = stringUrl.replace(AnnotationManager.CHR_LOCATION, region.start.chr.toNormalisedString());
			stringUrl = stringUrl.replace(AnnotationManager.START_LOCATION, region.start.bp.toString());
			stringUrl = stringUrl.replace(AnnotationManager.END_LOCATION, region.end.bp.toString());
			
			return stringUrl;
		} else {
			return "";
		}
		
	}
	
	protected void reportException(Exception e) {
		e.printStackTrace();
	}
	
	protected void showDialog(String title, String message, String details, boolean showDetails, boolean modal) {
		System.out.println("showDialog not implemented: " + title + "\t" +  message + "\t" + details + "\t" + modal);
	}
	
	protected void openExternalBrowser(String url) {
		System.out.println("openExternalBrowser not implemented: " + url);
	}
	
	protected void runBlockingTask(String taskName, Runnable runnable) {
		System.out.println("runBlockingTask not implemented: " + taskName + "\t" +  runnable);
	}
}
