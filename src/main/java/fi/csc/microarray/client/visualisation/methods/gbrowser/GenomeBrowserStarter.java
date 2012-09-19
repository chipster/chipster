package fi.csc.microarray.client.visualisation.methods.gbrowser;



/**
 * Quick and dirty starter utility for genome browser development and debugging.
 */
public class GenomeBrowserStarter {

//	private static final File BAM_DATA_FILE;
//	private static final File BAI_DATA_FILE;
//	private static final File CYTOBAND_FILE;
//	private static final File CYTOBAND_REGION_FILE;
//	private static final File GTF_ANNOTATION_FILE;
//
//	private static final String dataPath;
//	
//	static {
//		
//		dataPath = System.getProperty("user.home") + "/chipster/ohtu/";
//
//		BAM_DATA_FILE = new File(dataPath + "ohtu-within-chr.bam");
//		BAI_DATA_FILE = new File(dataPath + "ohtu-within-chr.bam.bai");
//		CYTOBAND_FILE = new File(dataPath + "Homo_sapiens.GRCh37.65.cytobands.txt");
//		CYTOBAND_REGION_FILE = new File(dataPath + "Homo_sapiens.GRCh37.65.seq_region.txt");
//		GTF_ANNOTATION_FILE = new File(dataPath + "Homo_sapiens.GRCh37.65.gtf");
//	}
//	
//	private static void checkData() {
//		File[] files = new File[] { BAM_DATA_FILE, BAI_DATA_FILE, CYTOBAND_FILE, CYTOBAND_REGION_FILE, GTF_ANNOTATION_FILE };
//		boolean fileNotFoundFail = false;
//		for (File file : files) {
//			if (!file.exists()) {
//				System.err.println("File not found: " + file);
//				fileNotFoundFail = true;
//			}
//		}
//		if (fileNotFoundFail) {
//			System.exit(1);
//		}
//	}
//
//	public static void main(String[] args) throws IOException {
//
//		JFrame frame = new JFrame();
//		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
//		getGenomeBrowserPanel();
//
//		frame.add(panel);
//		//frame.add(new JButton("Button"));
//		frame.setSize(PREVIEW_WIDTH, PREVIEW_HEIGHT);
//		
//		frame.setVisible(true);
//	}
//	
//	private static TooltipAugmentedChartPanel panel;
//	private static GenomePlot plot;
//	protected static int PREVIEW_WIDTH = 1280;
//	protected static int PREVIEW_HEIGHT = 768;
//	
//	public static TooltipAugmentedChartPanel getGenomeBrowserPanel() throws FileNotFoundException, MalformedURLException {
//		
//		checkData();
//		
//		boolean horizontal = true;
//		
//		panel = new TooltipAugmentedChartPanel();
//		
//		plot = new GenomePlot(panel, horizontal);
//		
//		TrackFactory.addCytobandTracks(plot, new CytobandDataSource(CYTOBAND_FILE, CYTOBAND_REGION_FILE));
//		
//		TrackFactory.addTitleTrack(plot, "Annotations");		
//		TrackFactory.addGeneTracks(plot, new LineDataSource(GTF_ANNOTATION_FILE));		
//		TrackFactory.addThickSeparatorTrack(plot);
//
//		TrackFactory.addReadTracks(
//				plot, 
//				new SAMDataSource(BAM_DATA_FILE, BAI_DATA_FILE),
//				SAMHandlerThread.class,
//				null,
//				"Control"
//		);
//		
//		TrackFactory.addReadTracks(
//				plot, 
//				new SAMDataSource(BAM_DATA_FILE, BAI_DATA_FILE),
//				SAMHandlerThread.class,
//				null,
//				"Treatment"
//		);
//		
//		panel.setChart(new JFreeChart(plot));
//		panel.setPreferredSize(new Dimension(PREVIEW_WIDTH, PREVIEW_HEIGHT));
//		panel.setCursor(new Cursor(Cursor.HAND_CURSOR));
//		
//		RegionDouble region = new RegionDouble(144151000d, 144154000d, new Chromosome("1"));
//		
//		plot.getDataView().setBpRegion(region, false);
//		
//		for (View view : plot.getViews()){
//			panel.addMouseListener(view);
//			panel.addMouseMotionListener(view);
//			panel.addMouseWheelListener(view);
//		}
//
//		return panel;
//	}
}
