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
import javax.swing.SwingUtilities;

import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.BedTabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.CytobandHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GeneSearchHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GtfTabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.IndexedFastaHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.SAMHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixSummaryHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.CytobandDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.IndexedFastaDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.SAMDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixSummaryDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationScrollGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserChartPanel;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserSettings;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GeneIndexActions;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.ScrollGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TooltipAugmentedChartPanel;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.ViewLimiter;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.AnnotationType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.GenomeAnnotation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.BedLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.ChromosomeBinarySearch;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.CnaConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.CnaLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.GtfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.GtfToFeatureConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.LineToRegionConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.RandomAccessLineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.VcfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackFactory;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.BrowserLauncher;
import fi.csc.microarray.util.IOUtils;

/**
 * Main class of genome browser visualisation. Depends on JFreeChart, SwingX, tribble, Picard and 
 * Chipster util package, but should not depend on any other Chipster code. All Chipster specific 
 * functionality is in class ChipsterGBrowserVisualisation.
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
		READS(true),
		HIDDEN(false), 
		VCF(true), 
		GTF(true),
		CNA_CALLS(true), 
		CNA_LOGRATIOS(true);

		public boolean isToggleable;

		private TrackType(boolean toggleable) {
			this.isToggleable = toggleable;
		}
	}
		
	public static class DataUrl {

		private URL url;
		private String name;

		public DataUrl(URL data, String name) {
			this.url = data;
			this.name = name;
		}

		public String getName() {
			return name;
		}

		public InputStream getInputStream() throws IOException, URISyntaxException {

			//Assume local
			return new FileInputStream(new File(url.toURI()));
		}

		public File getLocalFile() throws IOException, URISyntaxException {
			//Assume local
			return new File(url.toURI());
		}

		public URL getUrl() {
			return url;
		}
	}
	
	public static class Interpretation {
		
		private TrackType type;
		private List<DataUrl> summaryDatas = new LinkedList<DataUrl>();
		private DataUrl primaryData;
		private DataUrl indexData;
		private String name;

		public Interpretation(TrackType type, DataUrl primaryData) {
			this.type = type;
			this.primaryData = primaryData;
		}

		public TrackType getType() {
			return type;
		}

		public void setType(TrackType type) {
			this.type = type;
		}

		public List<DataUrl> getSummaryDatas() {
			return summaryDatas;
		}

		public void setSummaryDatas(List<DataUrl> summaryDatas) {
			this.summaryDatas = summaryDatas;
		}

		public DataUrl getPrimaryData() {
			return primaryData;
		}

		public void setPrimaryData(DataUrl primaryData) {
			this.primaryData = primaryData;
		}

		public DataUrl getIndexData() {
			return indexData;
		}

		public void setIndexData(DataUrl indexData) {
			this.indexData = indexData;
		}
		
		public void setName(String name) {
			this.name = name;
		}
		
		public String getName() {
			if (name != null) {
				return name;
			} else {
				return primaryData.getName();
			}
		}
	}

	public static class TrackDefinition {

		public Interpretation interpretation;
		public JCheckBox checkBox;
		public String name;
		public TrackGroup trackGroup = null;

		public TrackDefinition(String name, Interpretation interpretation) {
			this.name = name;
			this.interpretation = interpretation;
		}

		public void setTrackGroup(TrackGroup trackGroup) {
			this.trackGroup = trackGroup;
		}
	}
	
	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";
	
	private List<TrackDefinition> tracks = new LinkedList<TrackDefinition>();

	private GBrowserPlot plot;

	private JPanel plotPanel = new JPanel(new CardLayout());

	private AnnotationManager annotationManager;

	private GeneIndexActions gia;

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
		tracks.add(new TrackDefinition(AnnotationManager.AnnotationType.GTF_TABIX.getId(), new Interpretation(TrackType.GENES, null)));
		tracks.add(new TrackDefinition(AnnotationManager.AnnotationType.CYTOBANDS.getId(), new Interpretation(TrackType.CYTOBANDS, null)));


		for (int i = 0; i < interpretations.size(); i++) {
			Interpretation interpretation = interpretations.get(i);
			tracks.add(new TrackDefinition(interpretation.getName(), interpretation));
		}

		// update the dataset switches in the settings panel
		settings.updateDatasetSwitches();
	}

	public void setFullHeight(boolean fullHeight) {

		plot.setFullLayoutMode(fullHeight);		
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
	
	public void updateCoverageScale() {
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
	public void updateTracks() {

		//Remove tracks
		plot.getOverviewView().clean();
		plot.getDataView().clean();
		
		//There is a reference to track objects in scroll bars
		plot.chartPanel.clean();

		Genome genome = getGenome();
		
		ScrollGroup overview = new ScrollGroup("Overview");
		AnnotationScrollGroup annotations = new AnnotationScrollGroup();
		
		SeparatorTrack3D separator = new SeparatorTrack3D(0, Long.MAX_VALUE, true);
		separator.setView(plot.getDataView());
		plot.getDataView().addTrackGroup(new TrackGroup(separator));

		// Add selected annotation tracks
		for (TrackDefinition track : tracks) {
			if (track.checkBox.isSelected()) {
				switch (track.interpretation.type) {
				case CYTOBANDS:

					URL cytobandUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.CYTOBANDS);

					try {
						
						if (cytobandUrl != null) {
							CytobandDataSource cytobandDataSource = new CytobandDataSource(cytobandUrl);
							AreaRequestHandler cytobandRequestHandler = new CytobandHandlerThread(cytobandDataSource);

							overview.addTrackGroup(TrackFactory.getCytobandTrackGroup(plot, cytobandRequestHandler));

							this.viewLimiter = new ViewLimiter(plot.getOverviewView().getQueueManager(), 
									cytobandRequestHandler, plot.getOverviewView());
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

					URL gtfUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.GTF_TABIX);

					URL gtfIndexUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.GTF_TABIX_INDEX);

					URL repeatUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REPEAT);

					URL repeatIndexUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REPEAT_INDEX);

					AreaRequestHandler gtfRequestHandler = null;
					AreaRequestHandler repeatRequestHandler = null;
					
					try {
						if (gtfUrl != null && gtfIndexUrl != null) {
							TabixDataSource gtfDataSource = new TabixDataSource(gtfUrl, gtfIndexUrl);
							gtfRequestHandler = new GtfTabixHandlerThread(gtfDataSource);
							
							//Init gene search
							URL geneUrl = annotationManager.getAnnotation(
									genome, AnnotationManager.AnnotationType.GENE_CHRS).getUrl();
							LineDataSource geneDataSource = new LineDataSource(geneUrl, GeneSearchHandler.class);
							GeneSearchHandler geneRequestHandler = new GeneSearchHandler(geneDataSource);

							gia = new GeneIndexActions(plot.getDataView().getQueueManager(), gtfRequestHandler, geneRequestHandler);
							
						}

						if (repeatUrl != null && repeatIndexUrl != null) {
							TabixDataSource repeatDataSource = new TabixDataSource(repeatUrl, repeatIndexUrl);
							repeatRequestHandler = new BedTabixHandlerThread(repeatDataSource);
						}

						//Show ruler track even if there are now data sources
						TrackGroup geneGroup = TrackFactory.getGeneTrackGroup(plot, gtfRequestHandler, repeatRequestHandler, false);
						track.setTrackGroup(geneGroup);
						annotations.addTrackGroup(geneGroup);

					} catch (URISyntaxException e) {
						reportException(e);
					} catch (IOException e) {
						reportException(e);
					}
					break;

				case REFERENCE:
					// integrated into reads
					break;

				case TRANSCRIPTS:
					// integrated into genes
					break;
				default:
					break;
				}
			}
		}

		plot.getOverviewView().addScrollGroup(overview);		
		plot.getDataView().addScrollGroup(annotations);		
		plot.getDataView().addTrackGroup((TrackFactory.getThickSeparatorTrackGroup(plot)));
		ScrollGroup samples = new ScrollGroup("Samples", true);

		boolean firstReadTrack = true;
		
		// Add selected read tracks
		for (TrackDefinition track : tracks) {
			if (track.checkBox.isSelected()) {

				DataUrl dataUrl;
				try {
					dataUrl = track.interpretation.primaryData;
					AreaRequestHandler treatmentRequestHandler;
					if (track.interpretation.type == TrackType.READS) {
						
						if (!firstReadTrack) {
							samples.addTrackGroup((TrackFactory.getThinSeparatorTrackGroup(plot)));
						} else {
							firstReadTrack = false;
						}

						URL fastaUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REFERENCE);
						URL fastaIndexUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REFERENCE_INDEX);

						AreaRequestHandler refSeqRequestHandler = null;
						
						if (fastaUrl != null && fastaIndexUrl != null) {
							IndexedFastaDataSource refSeqDataSource = new IndexedFastaDataSource(fastaUrl, fastaIndexUrl);
							refSeqRequestHandler = new IndexedFastaHandlerThread(refSeqDataSource);
						}

						if (track.interpretation.summaryDatas.size() == 0) {
							// No precomputed summary data
							
							DataSource treatmentData = createReadDataSource(track.interpretation.primaryData, track.interpretation.indexData, tracks);
							treatmentRequestHandler = new SAMHandlerThread(treatmentData);

							TrackGroup readGroup = TrackFactory.getReadTrackGroup(
									plot, treatmentRequestHandler, 
									refSeqRequestHandler, 
									track.interpretation.primaryData.getName());

							track.setTrackGroup(readGroup);
							
							samples.addTrackGroup(readGroup);

						} else { 
							// Has precomputed summary data
							DataSource treatmentData = createReadDataSource(track.interpretation.primaryData, track.interpretation.indexData, tracks);
							treatmentRequestHandler = new SAMHandlerThread(treatmentData);
							
							DataSource symmaryData = new TabixDataSource(dataUrl.getUrl(), null);
							AreaRequestHandler summaryRequestHandler = new TabixSummaryHandlerThread(symmaryData);
							
							TrackGroup readGroupWithSummary = TrackFactory.getReadSummaryTrackGroup(
									plot, treatmentRequestHandler, refSeqRequestHandler, 
									track.interpretation.primaryData.getName(), summaryRequestHandler);
							track.setTrackGroup(readGroupWithSummary);
							samples.addTrackGroup(readGroupWithSummary);
						}
					}
				} catch (IOException e) {
					reportException(e);
				} catch (URISyntaxException e) {
					reportException(e);
				} catch (GBrowserException e) {
					reportException(e);
				}
			}
		}

		if (firstReadTrack) {// there wasn't any read tracks, add a separate reference  track

			//This track has fixed size now, layout system understands it only when the scrolling is disabled
			samples.setScrollEnabled(false);

			URL fastaUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REFERENCE);
			URL fastaIndexUrl = getAnnotationUrl(genome, AnnotationManager.AnnotationType.REFERENCE_INDEX);

			IndexedFastaDataSource refSeqDataSource = null;

			if (fastaUrl != null && fastaIndexUrl != null) {
				try {
					 
					refSeqDataSource = new IndexedFastaDataSource(fastaUrl, fastaIndexUrl);
					AreaRequestHandler refSeqRequestHandler = new IndexedFastaHandlerThread(refSeqDataSource);

					TrackGroup readGroup = TrackFactory.getReadTrackGroup(
							plot, null, 
							refSeqRequestHandler, 
							settings.getGenome().toString());

					samples.addTrackGroup(readGroup);

				} catch (URISyntaxException e) {
					e.printStackTrace();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

		plot.getDataView().addScrollGroup(samples);
		plot.getDataView().addTrackGroup(TrackFactory.getThickSeparatorTrackGroup(plot));
		ScrollGroup analysis = new ScrollGroup("Analysis", false);

		boolean firstPeakTrack = true;
		
		// Add selected peak tracks
		for (TrackDefinition track : tracks) {
			if (track.checkBox.isSelected()) {

				DataUrl dataUrl = track.interpretation.primaryData;
				
				//Add separators
				switch (track.interpretation.type) {
				case REGIONS:
				case VCF:
				case GTF:
				case CNA_CALLS:
				case CNA_LOGRATIOS:
					
					if (!firstPeakTrack) {
						analysis.addTrackGroup(TrackFactory.getThinSeparatorTrackGroup(plot));
					} else {
						firstPeakTrack = false;
					}
					break;
				default:
					break;
				}	
				
				switch (track.interpretation.type) {
				case REGIONS:
					
					analysis.addTrack(TrackFactory.getTitleTrack(plot, track.interpretation.primaryData.getName()));										

					try {						
						AreaRequestHandler conversion = new LineToRegionConversion(dataUrl.getUrl(), new BedLineParser(true));
						analysis.addTrackGroup(TrackFactory.getPeakTrackGroup(plot, conversion));
						
					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					}
					break;
					
				case VCF:

					analysis.addTrack(TrackFactory.getTitleTrack(plot, track.interpretation.primaryData.getName()));

					try {						
						AreaRequestHandler conversion = new LineToRegionConversion(dataUrl.getUrl(), new VcfLineParser());
						analysis.addTrackGroup(TrackFactory.getPeakTrackGroup(plot, conversion));
						
					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					}
					break;
				case GTF:

					analysis.addTrack(TrackFactory.getTitleTrack(plot, track.interpretation.primaryData.getName()));										
					analysis.setScrollEnabled(true);

					try {
						//DataSource gtfData = new LineDataSource(fileUrl, GtfToFeatureConversion.class);
						DataSource gtfData = new RandomAccessLineDataSource(dataUrl.getUrl());
						GtfToFeatureConversion gtfConversion = new GtfToFeatureConversion(gtfData, this);						
						analysis.addTrackGroup(TrackFactory.getGeneTrackGroup(plot, gtfConversion, null, true));
						
					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					} 
					break;
					
				case CNA_CALLS:
				case CNA_LOGRATIOS:

					analysis.addTrack(TrackFactory.getTitleTrack(plot, track.interpretation.primaryData.getName()));										
					analysis.setScrollEnabled(true);
					
					//Header has to be read to know the number of samples

					RandomAccessLineDataSource cnaData;

					try {
						cnaData = new RandomAccessLineDataSource(dataUrl.getUrl());
						CnaConversion conversion = new CnaConversion(cnaData, this);			
						
						cnaData.setLineReaderPosition(0);
						String header = cnaData.getNextLine();
						CnaLineParser parser = new CnaLineParser();
						parser.setLine(header);
						
						LinkedList<String> internalSampleNames = parser.getSampleNames();
						LinkedList<String> sampleNames = this.getSampleNames(internalSampleNames, dataUrl);
						
						boolean showCalls = (track.interpretation.type == TrackType.CNA_CALLS);
						boolean showLogratios = (track.interpretation.type == TrackType.CNA_LOGRATIOS);
						
						analysis.addTrackGroup(TrackFactory.getCnaTrackGroup(plot, conversion, sampleNames, showCalls, showLogratios));
						
					} catch (FileNotFoundException e) {
						reportException(e);
					} catch (URISyntaxException e) {
						reportException(e);
					} catch (IOException e) {
						reportException(e);
					} catch (GBrowserException e) {
						reportException(e);
					}						
					
					break;

				default:
					break;
				}				
			}
		}

		if (analysis.getTrackGroups().size() > 0) {
			plot.getDataView().addScrollGroup(analysis);
		}

		// End 3D effect
		SeparatorTrack3D separator2 = new SeparatorTrack3D(0, Long.MAX_VALUE, false);
		separator2.setView(plot.getDataView());
		plot.getDataView().addTrackGroup(new TrackGroup(separator2));
		
		//This does not fire area requests, but they are created separately when location is known, 
		//i.e. when the Go button is pressed or if dataset switches are used  
		plot.initializeTracks();
	}

	/**
	 * Create DataSource for SAM/BAM files
	 * 
	 * @param tracks
	 * 
	 * @param url
	 * @return
	 * @throws MicroarrayException
	 *             if index file is not selected properly
	 * @throws IOException
	 *             if opening data files fails
	 * @throws URISyntaxException 
	 * @throws GBrowserException 
	 */
	public DataSource createReadDataSource(DataUrl data, DataUrl indexData, List<TrackDefinition> tracks)
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

		}

		return dataSource;
	}

	public URL getAnnotationUrl(Genome genome, AnnotationManager.AnnotationType type) {
		GenomeAnnotation annotation = annotationManager.getAnnotation(
				genome, type);
		if (annotation != null) {
			return annotation.getUrl();					
		} else {
			return null;
		}
	}

	public void showVisualisation() {

		//Clean old data layers
		if (plot != null) {
			plot.clean();
		}

		// Create the chart panel with tooltip support				
		TooltipAugmentedChartPanel chartPanel = new TooltipAugmentedChartPanel();
		this.plot = new GBrowserPlot(chartPanel, true);
		
		((GBrowserChartPanel)chartPanel).setPlot(plot);

		//Set default location to plot to avoid trouble in track initialization. 
		plot.getDataView().setBpRegion(new RegionDouble(
				settings.getLocation() - settings.getViewSize() / 2.0, settings.getLocation() + settings.getViewSize() / 2.0, 
				settings.getChromosome()));
		
		plot.addDataRegionListener(settings);
				
		updateCoverageScale();
		
		updateTracks();
		
		settings.updateTracks();

		// Wrap GenomePlot in a panel
		chartPanel.setChart(new JFreeChart(plot));
		chartPanel.setCursor(new Cursor(Cursor.HAND_CURSOR));

		// Add mouse listeners
		for (GBrowserView view : plot.getViews()) {
			chartPanel.addMouseListener(view);
			chartPanel.addMouseMotionListener(view);
			chartPanel.addMouseWheelListener(view);
		}

		// Put panel on top of card layout
		if (plotPanel.getComponentCount() == 2) {
			plotPanel.remove(1);
		}

		setFullHeight(settings.isFullHeight());

		plotPanel.add(chartPanel, PLOTPANEL);
		plotPanel.addComponentListener(this);
		CardLayout cl = (CardLayout) (plotPanel.getLayout());
		cl.show(plotPanel, PLOTPANEL);
	}

	private GeneIndexActions getGeneIndexActions() {

		if (gia == null) {
			showDialog("Gene search failed", "Gene search is not initialized, is annotation data missing?", null, true, false, true, false);
		}
		return gia;
	}

	public void requestGeneSearch(String gene) {

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
									true, false, false, false);
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
					showDialog("Not found", "Gene was not found", null,	false, false, false, false);

				} else {

					// Update coordinate controls with gene's location

					Chromosome resultChr = new Chromosome(geneLocation.start.chr);

					if (settings.setChromosome(resultChr)) {

						setLocation(settings.getChromosome(), geneLocation.start.bp, geneLocation.end.bp);
					} else {
						showDialog("Different chromosome", 
								"Searched gene was found from chromosome " + resultChr + " but there is no data for that chromosome", "" + geneLocation, 
								true, false, false, false);
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
	
	public LinkedList<Chromosome> getChromosomeNames() throws IOException {

		// Gather all chromosome names from all indexed datasets (SAM/BAM)
		TreeSet<String> chromosomeNames = new TreeSet<String>(); 
		for (Interpretation interpretation : interpretations) {
			if (interpretation.type == TrackType.READS) {
				InputStream in = null;
				try {
					in  = interpretation.primaryData.getInputStream();
					chromosomeNames.addAll(SamBamUtils.readChromosomeNames(in));
				} catch (URISyntaxException e) {
					e.printStackTrace();
				} finally { 
					IOUtils.closeIfPossible(in);
				}
			}
		}

		// If we still don't have names, go through non-indexed datasets
		if (chromosomeNames.isEmpty()) {
			for (Interpretation interpretation : getInterpretations()) {
				
				boolean isBed = (interpretation.type == TrackType.REGIONS);
				boolean isVcf = (interpretation.type == TrackType.VCF);
				boolean isGtf = (interpretation.type == TrackType.GTF);
				boolean isCna = (interpretation.type == TrackType.CNA_CALLS);
				
				if (isBed || isVcf || isGtf || isCna) {
										
					try {
						
						DataUrl data = interpretation.primaryData;						
						ChromosomeBinarySearch chrSearch = null;
						
						if (isBed) {														
							chrSearch = new ChromosomeBinarySearch(data.getUrl(), new BedLineParser(true));														
						} else if (isVcf) {							
							chrSearch = new ChromosomeBinarySearch(data.getUrl(), new VcfLineParser());							
						} else if (isGtf) {
							chrSearch = new ChromosomeBinarySearch(data.getUrl(), new GtfLineParser());
						} else if (isCna) {
							chrSearch = new ChromosomeBinarySearch(data.getUrl(), new CnaLineParser());
						}
						
						for (Chromosome chr : chrSearch.getChromosomes()) {
							chromosomeNames.add(chr.toNormalisedString());
						}
						
					} catch (UnsortedDataException e) {
						this.showDialog("Unsorted data", e.getMessage(), null, true, false, true, true);					
						
					} catch (URISyntaxException e) {
						e.printStackTrace();
					} catch (GBrowserException e) {
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
	
	public String getExternalLinkUrl(AnnotationType browser) {
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
	
	public void openExternalBrowser(String url) {

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
	
	public List<TrackDefinition> getTracks() {
		return tracks;
	}
	
	/** 
	 * Override this method to customize error reporting
	 */
	public void reportException(Exception e) {
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
	public void showDialog(String title, String message, String details, boolean warning, boolean dialogShowDetails, boolean modal, boolean closeBrowser) {
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
	public void initialiseUserDatas() throws IOException {
		//Nothing to do if the files are already local
	}
	
	/** 
	 * Override this method to get the icons. Paths are defined in class GBrowserConstants.
	 */
	public ImageIcon getIcon(String path) {
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
	
	/**
	 * Override this convert internal sample names to prety names in phenodata
	 * 
	 * @param internalSampleNames
	 * @param dataUrl 
	 * @return
	 */
	public LinkedList<String> getSampleNames(
			LinkedList<String> internalSampleNames, DataUrl dataUrl) {
		return internalSampleNames;
	}

	/**
	 * Updates tracks to correspond with the settings and updates the data. This is used
	 * when the dataset visibility settings are changed and old location is shown with the new settings, whereas
	 * in initialization the tracks are created when the visualization opens, but data is requested only later after the "Go"
	 * button is pressed.
	 */
	public void updateData() {
		updateTracks();
		settings.updateVisibilityForTracks();
		plot.updateData();
	}
}
