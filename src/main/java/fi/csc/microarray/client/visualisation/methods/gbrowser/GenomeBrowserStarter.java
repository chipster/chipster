package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Cursor;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.visualisation.NonScalableChartPanel;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.SAMHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.CytobandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.GeneParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SequenceParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TranscriptParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.WIGParser;


/**
 * Quick started utility for developing genome browser. Might evolve into standalone version of
 * the browser some day.
 *
 */
public class GenomeBrowserStarter {

//	private static final File ELAND_DATA_FILE;
	private static final File BAM_DATA_FILE;
	private static final File BAI_DATA_FILE;

	private static final File MACS_DATA_FILE;
	private static final File WIG_DATA_FILE;
	private static final File URL_ROOT;
	private static final String annotationPath;
	
	static {
		
		annotationPath = "/home/" + System.getProperty("user.name") + "/chipster-share/";
		
//		ELAND_DATA_FILE = new File(annotationPath, "/ngs/STAT1/STAT1_treatment_aggregated_filtered_sorted_chr1.txt");
		BAM_DATA_FILE = new File(annotationPath + "/ngs/RNA-seq/pairedEnd_Berger/501Mel.sorted.bam");
		BAI_DATA_FILE = new File(annotationPath + "/ngs/RNA-seq/pairedEnd_Berger/501Mel.sorted.bam.bai");

		MACS_DATA_FILE = new File(annotationPath, "/ngs/STAT1/STAT1_peaks_sorted.bed");
		URL_ROOT = new File(annotationPath, "/ngs/annotations");
		
		WIG_DATA_FILE = new File(annotationPath, "/ngs/wig/GSM529979_chr1.wig.out");//variableStep - GSM545202.wig; fixedStep - GSM529979.wig
	}

	public static void main(String[] args) throws IOException {
		boolean horizontal = true;
        ChartPanel panel = new NonScalableChartPanel();
		GenomePlot plot = new GenomePlot(panel, horizontal);
		TrackFactory.addCytobandTracks(plot,
		        new ChunkDataSource(URL_ROOT, "Homo_sapiens.GRCh37.57_karyotype.tsv", new CytobandParser()));
		
		TrackFactory.addThickSeparatorTrack(plot);
		TrackFactory.addTitleTrack(plot, "Annotations");
		
		TrackFactory.addGeneTracks(plot,
		        new ChunkDataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_genes.tsv", new GeneParser()),
		        new ChunkDataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_transcripts.tsv", new TranscriptParser()),
		        new ChunkDataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_seq.tsv", new SequenceParser()));
//		TrackFactory.addMirnaTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_miRNA.tsv"));

		// Example peak: choromosome 21 in front of IFNAR2 gene (location 33,525,000)
		// Example peak: choromosome 21 in front of IFNAR1 gene (location 33,620,000)
		TrackFactory.addThickSeparatorTrack(plot);
		TrackFactory.addTitleTrack(plot, "WIG");
		
		TrackFactory.addWigTrack(plot,
		        new ChunkDataSource(WIG_DATA_FILE, new WIGParser()));
		
		TrackFactory.addThickSeparatorTrack(plot);
		TrackFactory.addTitleTrack(plot, "Peaks");
		
		TrackFactory.addPeakTrack(plot,
		        new ChunkDataSource(MACS_DATA_FILE, new BEDParser()));

		TrackFactory.addThickSeparatorTrack(plot);
		TrackFactory.addReadTracks(
				plot, 
				new SAMDataSource(BAM_DATA_FILE, BAI_DATA_FILE),
				SAMHandlerThread.class,
				new ChunkDataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_seq.tsv", new SequenceParser()),
				"Reads"
		);

		TrackFactory.addRulerTrack(plot);
		plot.start("1", 1024 * 1024 * 250d, 1000000L, 100000L);
		
		panel.setChart(new JFreeChart(plot));
		panel.setPreferredSize(new Dimension(800, 2000));
		panel.setCursor(new Cursor(Cursor.HAND_CURSOR));
		
		for (View view : plot.getViews()){
			panel.addMouseListener(view);
			panel.addMouseMotionListener(view);
			panel.addMouseWheelListener(view);
		}

		JFrame frame = new JFrame();
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
		frame.add(panel);
		frame.pack();
		frame.setVisible(true);
	}
	
}
