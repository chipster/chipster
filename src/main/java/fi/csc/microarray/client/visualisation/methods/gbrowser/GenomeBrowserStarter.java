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
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.CytobandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.GeneParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SequenceParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TranscriptParser;


/**
 * Quick started utility for developing genome browser. Might evolve into standalone version of
 * the browser some day.
 *
 */
public class GenomeBrowserStarter {

	private static final File ELAND_DATA_FILE = new File("/home/zukauska/chipster-share/ngs/STAT1/STAT1_treatment_aggregated_filtered_sorted_chr1.txt");
	private static final File MACS_DATA_FILE = new File("/home/zukauska/chipster-share/ngs/STAT1/STAT1_peaks_sorted.bed");
	private static final File URL_ROOT;

	static {
			URL_ROOT = new File("/home/zukauska/chipster-share/ngs/annotations");
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
		        new ChunkDataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_transcripts.tsv", new TranscriptParser()));
//		TrackFactory.addMirnaTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_miRNA.tsv"));

		// Example peak: choromosome 21 in front of IFNAR2 gene (location 33,525,000)
		// Example peak: choromosome 21 in front of IFNAR1 gene (location 33,620,000)
		
		TrackFactory.addThickSeparatorTrack(plot);
		TrackFactory.addTitleTrack(plot, "Peaks");
		
		TrackFactory.addPeakTrack(plot,
		        new ChunkDataSource(MACS_DATA_FILE, new BEDParser()));

		TrackFactory.addThickSeparatorTrack(plot);
		TrackFactory.addReadTracks(
				plot, 
				new ChunkDataSource(ELAND_DATA_FILE, new ElandParser()),
				ChunkTreeHandlerThread.class,
				new ChunkDataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_seq.tsv", new SequenceParser()),
				"Reads"
		);


		TrackFactory.addRulerTrack(plot);
		plot.start("1", 1024 * 1024 * 250d);
		plot.moveDataBpRegion(1000000L, 100000L);
		
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
