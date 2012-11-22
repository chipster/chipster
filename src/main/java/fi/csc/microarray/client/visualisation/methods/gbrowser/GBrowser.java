package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.CardLayout;
import java.awt.Cursor;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.BedTabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GeneSearchHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GtfTabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixSummaryHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.CytobandDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.IndexedFastaDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.SAMDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixSummaryDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParserWithCoordinateConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.HeaderTsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.VcfParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.AnnotationType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.GenomeAnnotation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.tools.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.tools.RegionOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.tools.SamBamUtils;
import fi.csc.microarray.client.visualisation.methods.gbrowser.tools.UnsortedDataException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackFactory;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.view.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.view.ViewLimiter;
import fi.csc.microarray.util.BrowserLauncher;
import fi.csc.microarray.util.IOUtils;

/**
 * Main class of genome browser visualisation. Depends on JFreeChart, SwingX, tribble, Picard and 
 * Chipster util package, but should not depend on any other Chipster code. All Chipster specific 
 * functionality must be in class ChipsterGBrowserVisualisation.
 * 
 * @author klemela
 */
public class GBrowser implements ComponentListener {
	
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

		public boolean isToggleable;

		private TrackType(boolean toggleable) {
			this.isToggleable = toggleable;
		}
	}
		
	public static class DataFile {

		private File file;

		public DataFile(File data) {
			this.file = data;
		}

		public String getName() {
			return file.getName();
		}

		public InputStream getInputStream() throws IOException {

			return new FileInputStream(file);
		}

		public File getLocalFile() throws IOException {
			//Assume local
			return file;
		}
	}
	
	public static class Interpretation {
		
		private TrackType type;
		private List<DataFile> summaryDatas = new LinkedList<DataFile>();
		private DataFile primaryData;
		private DataFile indexData;

		public Interpretation(TrackType type, DataFile primaryData) {
			this.type = type;
			this.primaryData = primaryData;
		}

		public TrackType getType() {
			return type;
		}

		public void setType(TrackType type) {
			this.type = type;
		}

		public List<DataFile> getSummaryDatas() {
			return summaryDatas;
		}

		public void setSummaryDatas(List<DataFile> summaryDatas) {
			this.summaryDatas = summaryDatas;
		}

		public DataFile getPrimaryData() {
			return primaryData;
		}

		public void setPrimaryData(DataFile primaryData) {
			this.primaryData = primaryData;
		}

		public DataFile getIndexData() {
			return indexData;
		}

		public void setIndexData(DataFile indexData) {
			this.indexData = indexData;
		}
	}

	public static class Track {

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
	
	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";
	
	private List<Track> tracks = new LinkedList<Track>();

	private GBrowserPlot plot;

	private JPanel plotPanel = new JPanel(new CardLayout());

	private AnnotationManager annotationManager;

	private GeneIndexActions gia;

	private JScrollPane verticalScroller;
	private ViewLimiter viewLimiter;
	protected boolean geneSearchDone;
	
	private GBrowserSettings settings;
	
	private List<Interpretation> interpretations;

	public void initialise() throws Exception {

		// initialize annotations
		this.annotationManager = new AnnotationManager(this);
		this.annotationManager.initialize();

		settings = new GBrowserSettings();
		settings.initialise(this);		
	}
	
	private void createAvailableTracks() {

		// for now just always add genes and cytobands
		tracks.add(new Track(AnnotationManager.AnnotationType.GTF_TABIX.getId(), new Interpretation(TrackType.GENES, null)));
		tracks.add(new Track(AnnotationManager.AnnotationType.CYTOBANDS.getId(), new Interpretation(TrackType.CYTOBANDS, null)));


		for (int i = 0; i < interpretations.size(); i++) {
			Interpretation interpretation = interpretations.get(i);
			tracks.add(new Track(interpretation.primaryData.getName(), interpretation));
		}

		// update the dataset switches in the settings panel
		settings.updateDatasetSwitches();
	}

	protected void setFullHeight(boolean fullHeight) {

		if (fullHeight) {
			verticalScroller.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
		} else {
			verticalScroller.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_NEVER);
		}

		plot.setFullHeight(fullHeight);		
	}

	public JComponent getVisualisation(List<Interpretation> interpretations) throws IOException {
		
		this.interpretations = interpretations;
		
		settings.updateInterpretations();

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
	
	protected void updateCoverageScale() {
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
	protected void updateTracks() {

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
				default:
					break;
				}
			}
		}

		// Add selected read tracks
		for (Track track : tracks) {
			if (track.checkBox.isSelected()) {

				File file;
				try {
					file = track.interpretation.primaryData == null ? null : track.interpretation.primaryData.getLocalFile();
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
					reportException(e);
				} catch (URISyntaxException e) {
					reportException(e);
				} catch (GBrowserException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
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
						file = track.interpretation.primaryData.getLocalFile();
						fileUrl = file.toURI().toURL();

					} catch (IOException e) {
						reportException(e);
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
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					} catch (IOException e) {
						reportException(e);
					} catch (UnsortedDataException e) {
						showDialog("Unsorted data", e.getMessage(), null, true, false, true);
					}
					break;
				case REGIONS_WITH_HEADER:
					TrackFactory.addThickSeparatorTrack(plot);
					TrackFactory.addTitleTrack(plot, track.interpretation.primaryData.getName());

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
					TrackFactory.addTitleTrack(plot, track.interpretation.primaryData.getName());

					try {
						regionData = new ChunkDataSource(fileUrl, new VcfParser(), ChunkTreeHandlerThread.class);
						TrackFactory.addPeakTrack(plot, regionData);

					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					}
					break;
				default:
					break;
				}
			}
		}

		// End 3D effect
		plot.getDataView().addTrack(new SeparatorTrack3D(plot.getDataView(), 0, Long.MAX_VALUE, false));
	}
	
	/**
	 * Create DataSource for SAM/BAM files
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
	 * @throws GBrowserException 
	 */
	public DataSource createReadDataSource(DataFile data, DataFile indexData, List<Track> tracks)
			throws IOException, URISyntaxException, GBrowserException {
		DataSource dataSource = null;

		// Convert data bean into file
		File file = data == null ? null : data.getLocalFile();

		URL fileUrl = file.toURI().toURL();

		if (data.getName().contains(".bam-summary")) {
			dataSource = new TabixSummaryDataSource(fileUrl);

		} else if (data.getName().contains(".bam") || data.getName().contains(".sam")) {
			File indexFile = indexData.getLocalFile();
			URL indexFileUrl = indexFile.toURI().toURL();
			dataSource = new SAMDataSource(fileUrl, indexFileUrl);

		} else {
			dataSource = new ChunkDataSource(fileUrl, new ElandParser(), ChunkTreeHandlerThread.class);
		}

		return dataSource;
	}

	protected URL getAnnotationUrl(Genome genome, AnnotationManager.AnnotationType type) {
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
		this.plot = new GBrowserPlot(chartPanel, true);
		plot.addDataRegionListener(settings);
		
		((GBrowserChartPanel)chartPanel).setGenomePlot(plot);

		//Set default location to plot to avoid trouble in track initialization. 
		plot.getDataView().setBpRegion(new RegionDouble(
				settings.getLocation() - settings.getViewSize() / 2.0, settings.getLocation() + settings.getViewSize() / 2.0, 
				settings.getChromosome()), true);
				
		updateCoverageScale();
		
		updateTracks();
		
		settings.updateTracks();

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

	protected void requestGeneSearch(String gene) {

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
									true, false, false);
						}

					});
				}
			}
		});

		getGeneIndexActions().requestLocation(gene, new GeneIndexActions.GeneLocationListener() {

			@Override
			public void geneLocation(Region geneLocation) {

				geneSearchDone = true;

				if (geneLocation == null) {

					// Move to last known location
					settings.processLocationPanelInput();

					// Tell the user 
					showDialog("Not found", "Gene was not found", null,	false, false, false);

				} else {

					// Update coordinate controls with gene's location

					Chromosome resultChr = new Chromosome(geneLocation.start.chr);

					if (settings.setChromosome(resultChr)) {

						setLocation(settings.getChromosome(), geneLocation.start.bp, geneLocation.end.bp);
					} else {
						showDialog("Different chromosome", 
								"Searched gene was found from chromosome " + resultChr + " but there is no data for that chromosome", "" + geneLocation, 
								true, false, false);
					}
				}
			}
		});
	}
	
	public void setLocation(Chromosome chr, Long start, Long end) {

		// Move to selected region
		
		settings.setChromosome(chr);

		if (end == null) {
			end = start;
		}
		
		settings.setCoordinateFields((end + start) / 2, (end - start) * 2);

		// Update
		plot.moveDataBpRegion((Chromosome) settings.getChromosome(),
				settings.getLocation(), settings.getViewSize());

		// Set scale of profile track containing reads information
		this.plot.setReadScale(settings.getCoverageScale());
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
	
	protected LinkedList<Chromosome> getChromosomeNames() throws IOException {

		// Gather all chromosome names from all indexed datasets (SAM/BAM)
		TreeSet<String> chromosomeNames = new TreeSet<String>(); 
		for (Interpretation interpretation : interpretations) {
			if (interpretation.type == TrackType.READS) {
				InputStream in = null;
				try {
					in  = interpretation.primaryData.getInputStream();
					chromosomeNames.addAll(SamBamUtils.readChromosomeNames(in));
				} finally { 
					IOUtils.closeIfPossible(in);
				}
			}
		}

		// If we still don't have names, go through non-indexed datasets
		if (chromosomeNames.isEmpty()) {
			for (Interpretation interpretation : getInterpretations()) {
				if (interpretation.type == TrackType.REGIONS) {
					DataFile data = interpretation.primaryData;
					File file = data.getLocalFile();
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

		return chromosomes;
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
		
		//FIXME remove if works without this
//		//Move to last location
//		settings.processLocationPanelInput();
		plot.redraw();
	}

	public List<Interpretation> getInterpretations() {
		return interpretations;
	}
	
	protected String getExternalLinkUrl(AnnotationType browser) {
		settings.getGenome();
		URL url = annotationManager.getAnnotation(settings.getGenome(), browser).getUrl();

		if (url != null && plot != null && plot.getDataView() != null && plot.getDataView().getBpRegion() != null) {
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
	
	protected void openExternalBrowser(String url) {

		try {
			BrowserLauncher.openURL(url);
		} catch (Exception e) {
			reportException(e);
		}
	}
	
	public JPanel getParameterPanel() {
		return settings.getParameterPanel();
	}
	
	public AnnotationManager getAnnotationManager() {
		return annotationManager;
	}
	
	public List<Track> getTracks() {
		return tracks;
	}
	
	/** 
	 * Override this method to customize error reporting
	 */
	protected void reportException(Exception e) {
		e.printStackTrace();
	}
	

	/**
	 * Override this method to show custom dialogs
	 * 
	 * @param title
	 * @param message
	 * @param details
	 * @param warning true for "warning" level message, false for "info"
	 * @param dialogShowDetails Show details by default
	 * @param modal
	 */
	protected void showDialog(String title, String message, String details, boolean warning, boolean dialogShowDetails, boolean modal) {
		System.out.println("showDialog not implemented: " + title + "\t" +  message + "\t" + details);
	}
	
	/** 
	 * Override this method to lock the gui during heavy tasks
	 */
	public void runBlockingTask(String taskName, Runnable runnable) {
		System.out.println("runBlockingTask: " + taskName);
		new Thread(runnable).start();
	}

	/** 
	 * Override this method to perform any slow initialization of the data files
	 */
	protected void initialiseUserDatas() throws IOException {
		//Nothing to do if the files are already local
	}
	
	/** 
	 * Override this method to get the icons. Paths are defined in class GBrowserConstants.
	 */
	protected ImageIcon getIcon(String path) {
		System.out.println("getIcon not implemented");
		return new ImageIcon();
	}

	/** 
	 * Override this method to let user decide if s/he want's to download to annotations.
	 */
	public void openDownloadAnnotationsDialog(Genome genome) {
		//Don't ask, just do it
		try {
			getAnnotationManager().downloadAnnotations(genome);
		} catch (IOException e) {
			reportException(e);
		}
	}
	
	/**
	 * Override this method to specify location for remote annotations
	 */
	public URL getRemoteAnnotationsUrl() throws Exception {
		//"http://chipster-filebroker.csc.fi:8080/public/annotations/"
		System.out.println("getRemoteAnnotationsUrl not implemented");
		return null;
	}
	
	/** 
	 * Override this method to specify location for local annotations.
	 * 
	 * WARNING: FILES IN THIS FOLDER WILL BE DELETED!
	 */
	public File getLocalAnnotationDir() throws IOException {
		//"~/.chipster/annotations/"
		System.out.println("getLocalAnnotationDir not implemented");
		return null;
	}
}
