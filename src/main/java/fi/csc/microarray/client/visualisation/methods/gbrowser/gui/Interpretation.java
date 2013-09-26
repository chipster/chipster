package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.BamDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.BamToCoverageConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.BamToCoverageEstimateConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.BamToDetailsConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.GtfToFeatureConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.IndexedFastaConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.BedLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.BedTabixToRegionConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ChromosomeBinarySearch;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.CnaConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.CnaLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.CytobandConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.FileLineConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.GeneSearchConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.GtfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.RandomAccessLineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ScatterplotFileLineConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.TsvLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.VcfLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;


public class Interpretation {
	
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
		CNA(true), 		
		TSV(true),
		TSV_WITH_ROW_ID(true);

		public boolean isToggleable;

		private TrackType(boolean toggleable) {
			this.isToggleable = toggleable;
		}
	}
	
	private TrackType type;
	private DataUrl primaryData;
	private DataUrl indexData;
	private String name;
	private CnaConversion cnaDataThread;	
	private GtfToFeatureConversion gtfDataThread;
	private FileLineConversion vcfDataThread;
	private ScatterplotFileLineConversion bedLineDataThread;
	private FileLineConversion tsvDataThread;
	private TreeSet<Chromosome> chrNames;

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
	
	public static CytobandConversion getCytobandDataThread(GBrowser browser) {

		DataUrl cytobandUrl = browser.getAnnotationUrl(browser.getGenome(), AnnotationManager.AnnotationType.CYTOBANDS);

		if (cytobandUrl != null) {

			return new CytobandConversion(cytobandUrl, browser);
		}
		throw new IllegalStateException("Cytoband url is null");
	}
	
	public static GeneIndexActions getGeneSearchDataThread(GBrowser browser) {

		Genome genome = browser.getGenome();

		GtfToFeatureConversion gtfDataThread = getAnnotationDataThread(browser);

		if (gtfDataThread != null) {
			//Init gene search
			DataUrl geneUrl = browser.getAnnotationManager().getAnnotation(
					genome, AnnotationManager.AnnotationType.GENE_CHRS).getUrl();

			GeneSearchConversion geneRequestHandler = new GeneSearchConversion(geneUrl, browser);

			return new GeneIndexActions(browser.getPlot().getDataView().getQueueManager(), gtfDataThread, geneRequestHandler);
		} 
		
		throw new IllegalStateException("Can't initialize gene search without gtf data");
	}	
	
	public static GtfToFeatureConversion getAnnotationDataThread(GBrowser browser) {
		
		Genome genome = browser.getGenome();
		
		DataUrl gtfUrl = browser.getAnnotationUrl(genome, AnnotationManager.AnnotationType.GTF);
		DataUrl gtfTabixUrl = browser.getAnnotationUrl(genome, AnnotationManager.AnnotationType.GTF_TABIX);
		DataUrl gtfIndexUrl = browser.getAnnotationUrl(genome, AnnotationManager.AnnotationType.GTF_TABIX_INDEX);
		
		GtfToFeatureConversion gtfDataThread = null;

		if (gtfUrl != null) {			
			gtfDataThread = new GtfToFeatureConversion(gtfUrl, null, browser);
			
		} else 	if (gtfTabixUrl != null && gtfIndexUrl != null) {
			gtfDataThread = new GtfToFeatureConversion(gtfTabixUrl, gtfIndexUrl, browser);
		}

		return gtfDataThread;
	}
	
	public static BedTabixToRegionConversion getRepeatDataThread(GBrowser browser) {
		
		Genome genome = browser.getGenome();

		DataUrl repeatUrl = browser.getAnnotationUrl(genome, AnnotationManager.AnnotationType.REPEAT);
		DataUrl repeatIndexUrl = browser.getAnnotationUrl(genome, AnnotationManager.AnnotationType.REPEAT_INDEX);
					
		if (repeatUrl != null && repeatIndexUrl != null) {							
			return new BedTabixToRegionConversion(repeatUrl, repeatIndexUrl, browser);
		}
		
		return null;
	}
	
	public static IndexedFastaConversion getReferenceDataThread(GBrowser browser) {
		
		Genome genome = browser.getGenome();

		DataUrl fastaUrl = browser.getAnnotationUrl(genome, AnnotationManager.AnnotationType.REFERENCE);
		DataUrl fastaIndexUrl = browser.getAnnotationUrl(genome, AnnotationManager.AnnotationType.REFERENCE_INDEX);

		IndexedFastaConversion refSeqDataThread = null;
		
		if (fastaUrl != null && fastaIndexUrl != null) {
			refSeqDataThread = new IndexedFastaConversion(fastaUrl, fastaIndexUrl, browser);
		}
		
		return refSeqDataThread;
	}
	

	public BamToDetailsConversion getBamDetailsDataThread(GBrowser browser) {

		if (getType() == TrackType.READS) {

			BamDataSource dataSource;
			try {
				//Create always a new data source, because picard doesn't support concurrent access
				dataSource = new BamDataSource(getPrimaryData(), getIndexData());
				return new BamToDetailsConversion(dataSource, browser);
				
			} catch (URISyntaxException | IOException e) {
				browser.reportException(e);
			}
		}
		throw new IllegalStateException("requested DataThread is not compatible with the Interpreation type: " + getType());
	}

	public BamToCoverageConversion getBamCoverageDataThread(GBrowser browser) {

		if (getType() == TrackType.READS) {

			try {
				//Create always a new data source, because picard doesn't support concurrent access
				BamDataSource dataSource;
				dataSource = new BamDataSource(getPrimaryData(), getIndexData());
				return new BamToCoverageConversion(dataSource, browser);
				
			} catch (URISyntaxException | IOException e) {
				browser.reportException(e);
			}
		}
		throw new IllegalStateException("requested DataThread is not compatible with the Interpreation type: " + getType());
	}

	public BamToCoverageEstimateConversion getBamCoverageEstimateDataThread(GBrowser browser) {

		if (getType() == TrackType.READS) {
			
			BamDataSource dataSource;
			try {
				//Create always a new data source, because picard doesn't support concurrent access
				dataSource = new BamDataSource(getPrimaryData(), getIndexData());
				return new BamToCoverageEstimateConversion(dataSource, browser);
				
			} catch (URISyntaxException | IOException e) {
				browser.reportException(e);
			}
		}
		throw new IllegalStateException("requested DataThread is not compatible with the Interpreation type: " + getType());
	}
	
	public ScatterplotFileLineConversion getBedLineDataThread(GBrowser browser) {

		if (getType() == TrackType.REGIONS) {

			if (bedLineDataThread == null) {
				try {
					bedLineDataThread = new ScatterplotFileLineConversion(getPrimaryData(), new BedLineParser(true), browser);

				} catch (URISyntaxException | IOException e) {
					browser.reportException(e);
				} catch (GBrowserException e) {
					e.printStackTrace();
				}
			}
			return bedLineDataThread;
		}
		throw new IllegalStateException("requested DataThread is not compatible with the Interpreation type: " + getType());
	}
	
	public FileLineConversion getTsvDataThread(GBrowser browser) {

		if (getType() == TrackType.TSV || getType() == TrackType.TSV_WITH_ROW_ID) {

			if (tsvDataThread == null) {
				try {
					tsvDataThread = new FileLineConversion(getPrimaryData(), new TsvLineParser(getPrimaryData(), getType()), browser);

				} catch (URISyntaxException | IOException e) {
					browser.reportException(e);
				}
			}
		}
		return tsvDataThread;
	}
	
	public FileLineConversion getVcfDataThread(GBrowser browser) {

		if (getType() == TrackType.VCF) {

			if (vcfDataThread == null) {
				try {
					vcfDataThread = new FileLineConversion(getPrimaryData(), new VcfLineParser(), browser);

				} catch (URISyntaxException | IOException e) {
					browser.reportException(e);
				}
			}
			return vcfDataThread;
		}
		throw new IllegalStateException("requested DataThread is not compatible with the Interpreation type: " + getType());
	}
	
	public GtfToFeatureConversion getGtfDataThread(GBrowser browser) {

		if (getType() == TrackType.GTF) {

			if (gtfDataThread == null) {
				gtfDataThread = new GtfToFeatureConversion(getPrimaryData(), null, browser);
			}
			return gtfDataThread;
		}
		
		throw new IllegalStateException("requested DataThread is not compatible with the Interpreation type: " + getType());
	}
	
	public CnaConversion getCnaDataThread(GBrowser browser) {

		if (getType() == TrackType.CNA) {
			
			if (cnaDataThread == null) {

				try {
					cnaDataThread = new CnaConversion(new RandomAccessLineDataSource(getPrimaryData()), browser);
					
				} catch (URISyntaxException | IOException e) {
					browser.reportException(e);
				}									
			}
			return cnaDataThread;
		}
		throw new IllegalStateException("requested DataThread is not compatible with the Interpreation type: " + getType());
	}

	public TreeSet<Chromosome> getChromosomeNames() throws URISyntaxException, IOException, UnsortedDataException, GBrowserException {
		
		TreeSet<Chromosome> chromosomes = new TreeSet<>();
		
		if (getType() == TrackType.READS) {

			URL bam  = getPrimaryData().getUrl();
			URL index  = getIndexData().getUrl();

			for (String string : SamBamUtils.readChromosomeNames(bam, index)) {

				chromosomes.add(new Chromosome(string));
			}
		}		
		
		boolean isBed = (getType() == TrackType.REGIONS);
		boolean isVcf = (getType() == TrackType.VCF);
		boolean isGtf = (getType() == TrackType.GTF);
		boolean isTsv = (getType() == TrackType.TSV || getType() == TrackType.TSV_WITH_ROW_ID);
		boolean isCna = (getType() == TrackType.CNA);

		if (isBed || isVcf || isGtf || isCna || isTsv) {

			DataUrl data = getPrimaryData();						
			ChromosomeBinarySearch chrSearch = null;

			if (isBed) {					
				chrSearch = new ChromosomeBinarySearch(data, new BedLineParser(true));			
			} else if (isVcf) {							
				chrSearch = new ChromosomeBinarySearch(data, new VcfLineParser());				
			} else if (isGtf) {
				chrSearch = new ChromosomeBinarySearch(data, new GtfLineParser());				
			} else if (isCna) {
				chrSearch = new ChromosomeBinarySearch(data, new CnaLineParser());
			} else if (isTsv) {
				chrSearch = new ChromosomeBinarySearch(data, new TsvLineParser(data, getType()));
			}

			chrNames = chrSearch.getChromosomes();
			chromosomes.addAll(chrNames);
		}
		
		return chromosomes;
	}
}