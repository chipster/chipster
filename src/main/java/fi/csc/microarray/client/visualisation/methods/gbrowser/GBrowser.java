package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.CardLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.BamToCoverageConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.BamToCoverageEstimateConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.BamToDetailsConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.GtfToFeatureConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager.AnnotationType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager.GenomeAnnotation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserSettings;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GeneIndexActions;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Interpretation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Interpretation.TrackType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.ScrollGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.SelectionManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.ViewLimiter;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.CnaConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ScatterplotFileLineConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.AnnotationTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.BedTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CnaCallsTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CnaFrequenciesTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CnaLogRatiosTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CytobandTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SampleTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.SeparatorTrack3D;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackFactory;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.RegionTrackGroup;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.util.BrowserLauncher;

/**
 * Main class of genome browser visualisation. Depends on SwingX, tribble, Picard and 
 * Chipster util package, but should not depend on any other Chipster code. All Chipster specific 
 * functionality is in class ChipsterGBrowserVisualisation.
 * 
 * @author klemela
 */
public class GBrowser {

	final static String WAITPANEL = "waitpanel";
	final static String PLOTPANEL = "plotpanel";

	private GBrowserPlot plot;

	private JPanel plotPanel = new JPanel(new CardLayout());

	private AnnotationManager annotationManager;

	private GeneIndexActions gia;

	private ViewLimiter viewLimiter;
	protected boolean geneSearchDone;

	private GBrowserSettings settings;

	private List<Interpretation> interpretations;
	private LinkedList<String> sampleNames;
	private SelectionManager selectionManager;

	public void initialise() throws Exception {

		// initialize annotations
		this.annotationManager = new AnnotationManager(this);
		this.annotationManager.initialize();

		settings = new GBrowserSettings();
		settings.initialise(this);
		
		this.selectionManager = new SelectionManager(this);
	}

	public JComponent getVisualisation(List<Interpretation> interpretations) throws IOException {

		this.interpretations = interpretations;

		settings.updateInterpretations();				

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

	public Genome getGenome() {
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

		GBrowserView dataView = plot.getDataView();
		GBrowserView overviewView = plot.getOverviewView();
		
		//Remove tracks		
		overviewView.clean();
		dataView.clean();

		ScrollGroup overviewScrollGroup = new ScrollGroup("Overview");		
		createOverviewTracks(dataView, overviewView, overviewScrollGroup);	
		overviewView.addScrollGroup(overviewScrollGroup);
		
		SeparatorTrack3D separator = new SeparatorTrack3D(true);		
		separator.setView(dataView);
		dataView.addTrackGroup(new TrackGroup(separator));
		
		ScrollGroup annotations = new ScrollGroup("Annotations", true);
		createAnnotationTracks(dataView, annotations);			
		dataView.addScrollGroup(annotations);		
		
		dataView.addTrackGroup((TrackFactory.getThickSeparatorTrackGroup(plot)));
		
		ScrollGroup samples = new ScrollGroup("Samples", true);
		createSampleTracks(dataView, samples);
		dataView.addScrollGroup(samples);

		ScrollGroup analyses = new ScrollGroup("Analyses", false);
		createAnalysisTracks(dataView, analyses);			
		dataView.addScrollGroup(analyses);

		// End 3D effect
		SeparatorTrack3D separator2 = new SeparatorTrack3D(false);
		separator2.setView(dataView);
		dataView.addTrackGroup(new TrackGroup(separator2));

		//This does not fire data requests, but they are created separately when location is known, 
		//i.e. when the Go button is pressed or if dataset switches are used  
		plot.initializeDataResultListeners();
	}

	private void createAnalysisTracks(GBrowserView dataView,
			ScrollGroup analyses) {
		
		boolean firstAnalysisTrack = true;

		//Add separators between analysis tracks
		for (Interpretation interpretation : interpretations) {
						
			switch (interpretation.getType()) {
			case REGIONS:
			case VCF:
			case GTF:
			case CNA:
				
				if (!firstAnalysisTrack) {
					analyses.addTrackGroup(TrackFactory.getThinSeparatorTrackGroup(plot));
				} else {
					firstAnalysisTrack = false;
				}
				break;
				
			default:
				break;
			}				
			
			String title = getTitle(interpretation);			

			// Add selected analysis tracks
			switch (interpretation.getType()) {
			case REGIONS:			
				
				ScatterplotFileLineConversion bed = interpretation.getBedLineDataThread(this);
				TrackGroup bedTrackGroup = new BedTrackGroup(dataView, bed, title);				
				analyses.addTrackGroup(bedTrackGroup);

				break;

			case VCF:
				DataThread vcf = interpretation.getVcfDataThread(this);
				TrackGroup vcfTrackGroup = new RegionTrackGroup(dataView, vcf, title);
				analyses.addTrackGroup(vcfTrackGroup);
				break;
			case TSV:
				DataThread tsv = interpretation.getTsvDataThread(this);
				analyses.addTrackGroup(new RegionTrackGroup(dataView, tsv, title));
				break;				
			case TSV_WITH_ROW_ID:
				tsv = interpretation.getTsvDataThread(this);
				analyses.addTrackGroup(new RegionTrackGroup(dataView, tsv, title));
				break;
			case GTF:

				GtfToFeatureConversion gtfConversion = interpretation.getGtfDataThread(this);
				TrackGroup geneGroup = new AnnotationTrackGroup(dataView, gtfConversion, null, true);
				
				analyses.setScrollEnabled(true);
				geneGroup.setSettingsEnabled(true);
				
				analyses.addTrackGroup(geneGroup);
				
				break;

			case CNA:

				analyses.setScrollEnabled(true);
				
				CnaConversion cnaData = interpretation.getCnaDataThread(this);

				LinkedList<String> internalSampleNames = cnaData.getSampleNames();
				this.sampleNames = this.getSampleNames(internalSampleNames, interpretation.getPrimaryData());

				CnaFrequenciesTrackGroup frequencies = new CnaFrequenciesTrackGroup(dataView, cnaData, sampleNames, title);
				CnaCallsTrackGroup calls = new CnaCallsTrackGroup(dataView, cnaData, sampleNames, title);
				CnaLogRatiosTrackGroup logRatios = new CnaLogRatiosTrackGroup(dataView, cnaData, sampleNames, title);

				analyses.addTrackGroup(frequencies);
				analyses.addTrackGroup(TrackFactory.getThinSeparatorTrackGroup(plot));
				analyses.addTrackGroup(calls);
				analyses.addTrackGroup(TrackFactory.getThinSeparatorTrackGroup(plot));
				analyses.addTrackGroup(logRatios);
				break;

			default:
				break;
			}				
		}
		
		if (analyses.getTrackGroups().size() > 0) {
			dataView.addTrackGroup(TrackFactory.getThickSeparatorTrackGroup(plot));			
		}
	}

	private void createSampleTracks(GBrowserView dataView, ScrollGroup samples) {
		boolean firstReadTrack = true;

		// Add selected read tracks
		for (Interpretation interpretation : interpretations) {
			
			if (interpretation.getType() == TrackType.READS) {
				String title = getTitle(interpretation);

				if (!firstReadTrack) {
					samples.addTrackGroup((TrackFactory.getThinSeparatorTrackGroup(plot)));
				} else {
					firstReadTrack = false;
				}

				//A separate data thread for each sample to avoid concurrency
				DataThread refSeqRequestHandler = Interpretation.getReferenceDataThread(this);											

				BamToDetailsConversion details = interpretation.getBamDetailsDataThread(this);
				BamToCoverageConversion coverage = interpretation.getBamCoverageDataThread(this);
				BamToCoverageEstimateConversion estimate = interpretation.getBamCoverageEstimateDataThread(this);
				
				SampleTrackGroup readGroup = new SampleTrackGroup(dataView, details, coverage, estimate, refSeqRequestHandler, title);
				readGroup.initialise();

				samples.addTrackGroup(readGroup);
			}
		}

		if (firstReadTrack) {// there wasn't any read tracks, add a separate reference  track

			//This track has fixed size now
			samples.setScrollEnabled(false);

			DataThread refSeqRequestHandler = Interpretation.getReferenceDataThread(this);

			if (refSeqRequestHandler != null) {
				
				SampleTrackGroup readGroup = new SampleTrackGroup(
						dataView, null, null, null, refSeqRequestHandler, settings.getGenome().toString());
				
				readGroup.initialise();

				samples.addTrackGroup(readGroup);
			}
		}
	}

	private void createAnnotationTracks(GBrowserView dataView,
			ScrollGroup annotations) {
		
		DataThread gtfRequestHandler = Interpretation.getAnnotationDataThread(this);
		DataThread repeatRequestHandler = Interpretation.getRepeatDataThread(this);

		gia = Interpretation.getGeneSearchDataThread(this);

		TrackGroup geneGroup = new AnnotationTrackGroup(dataView, gtfRequestHandler, repeatRequestHandler, false);
		annotations.addTrackGroup(geneGroup);
	}

	private void createOverviewTracks(GBrowserView dataView,
			GBrowserView overviewView, ScrollGroup overviewScrollGroup) {
		
		DataThread cytobandDataThread = Interpretation.getCytobandDataThread(this);

		if (cytobandDataThread != null) {
								
			CytobandTrack overviewCytobands = new CytobandTrack(false);
			overviewCytobands.setView(overviewView);
			overviewCytobands.addDataThread(cytobandDataThread);
			
			TrackGroup cytobandTrackGroup =  new TrackGroup(overviewCytobands);
			overviewScrollGroup.addTrackGroup(cytobandTrackGroup);

			//View limits are based on cytoband data
			this.viewLimiter = new ViewLimiter(overviewView.getQueueManager(), 
					cytobandDataThread, overviewView);
			dataView.setViewLimiter(viewLimiter);
			overviewView.setViewLimiter(viewLimiter);
		}
	}

	private String getTitle(Interpretation interpretation) {
		if (interpretation != null && interpretation.getPrimaryData() != null) {
			return interpretation.getPrimaryData().getName();
		}
		return null;
	}

	public DataUrl getAnnotationUrl(Genome genome, AnnotationManager.AnnotationType type) {
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

		Region defaultLocation = new Region(
				(long)(settings.getLocation() - settings.getViewSize() / 2.0), 
				(long)(settings.getLocation() + settings.getViewSize() / 2.0), 
				settings.getChromosome());
				
		this.plot = new GBrowserPlot(this, true, defaultLocation);

		plot.addDataRegionListener(settings);

		updateCoverageScale();

		updateTracks();

		settings.updateTracks();

		// Put panel on top of card layout
		if (plotPanel.getComponentCount() == 2) {
			plotPanel.remove(1);
		}

		plotPanel.add(plot.getComponent(), PLOTPANEL);

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

		plotPanel.removeAll();

		if (plot != null) {
			plot.clean();
			plot = null;
		}

		//Remove references to tracks and data to free memory, even if the (hidden) parameter panel keeps actionListener
		//references to this object preventing garbage collection (when visualization is changed to none)
		if (interpretations != null) {
			interpretations.clear();
		}
		gia = null;	
	}

	public LinkedList<Chromosome> getChromosomeNames() throws IOException {

		// Gather all chromosome names from all indexed datasets (SAM/BAM)
		TreeSet<Chromosome> chromosomes = new TreeSet<>(); 
		try {
			for (Interpretation interpretation : interpretations) {
				if (interpretation.getType() == TrackType.READS) {

					chromosomes.addAll(interpretation.getChromosomeNames());
				}
			}

			// If we still don't have names, go through non-indexed datasets
			if (chromosomes.isEmpty()) {
				for (Interpretation interpretation : getInterpretations()) {
					if (interpretation.getType() != TrackType.READS) {	
						chromosomes.addAll(interpretation.getChromosomeNames());
					}
				}
			}
		} catch (URISyntaxException	| GBrowserException e) {
			reportException(e);
		}

		LinkedList<Chromosome> list = new LinkedList<Chromosome>();

		for (Chromosome chromosome : chromosomes) {
			list.add(chromosome);
		}

		return list;
	}

	public List<Interpretation> getInterpretations() {
		return interpretations;
	}

	public String getExternalLinkUrl(AnnotationType browser) {
		settings.getGenome();
		GenomeAnnotation urlAnnotation = annotationManager.getAnnotation(settings.getGenome(), browser);

		if (urlAnnotation != null) {
			URL url = null;
			try {
				url = urlAnnotation.getUrl().getUrl();
			} catch (IOException e) {
				//Just disable link
			}

			if (url != null && plot != null && plot.getDataView() != null && plot.getDataView().getBpRegion() != null) {
				String stringUrl = url.toString();
				Region region = plot.getDataView().getBpRegion();
				stringUrl = stringUrl.replace(AnnotationManager.CHR_LOCATION, region.start.chr.toNormalisedString());
				stringUrl = stringUrl.replace(AnnotationManager.START_LOCATION, region.start.bp.toString());
				stringUrl = stringUrl.replace(AnnotationManager.END_LOCATION, region.end.bp.toString());

				return stringUrl;
			} 
		}
		return "";
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
	@Deprecated
	public URL getRemoteAnnotationsUrl() throws Exception {
		//"http://chipster-filebroker.csc.fi:8080/public/annotations/"
		System.out.println("getRemoteAnnotationsUrl not implemented");
		return null;
	}

	/**
	 * Override this method to specify location for remote annotations
	 */
	public List<URL> getRemoteAnnotationFiles() throws Exception {

		System.out.println("getRemoteAnnotationFiles not implemented");
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
	 * Override this convert internal sample names to pretty names of the phenodata
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

	public GBrowserPlot getPlot() {
		return plot;
	}

	public SelectionManager getSelectionManager() {
		return selectionManager;
	}

	public void initializeDataResultListeners() {
		plot.initializeDataResultListeners();
		//This happens in initialization		
		if (gia != null) {
			gia.initializeDataResultListeners();
		}
	}
}
