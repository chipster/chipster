package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Cursor;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;



public class GenomeBrowserStarter {

	private static final File ELAND_DATA_FILE = new File("/home/klemela/chipster-share/ngs/STAT1/STAT1_treatment_aggregated_filtered_sorted_chr1.txt");
	private static final File MACS_DATA_FILE = new File("/home/klemela/chipster-share/ngs/STAT1/STAT1_peaks_sorted.bed");
	private static final File URL_ROOT;

	static {
			URL_ROOT = new File("/home/klemela/chipster-share/ngs/annotations");
	}

	public static void main(String[] args) throws IOException {
		GenomePlot plot = new GenomePlot(true);
		TrackFactory.addCytobandTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.57_karyotype.tsv"));
		TrackFactory.addGeneTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_genes.tsv"), new DataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_transcripts.tsv"));
//		TrackFactory.addMirnaTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_miRNA.tsv"));

		// Example peak: choromosome 21 in front of IFNAR2 gene (location 33,525,000)
		// Example peak: choromosome 21 in front of IFNAR1 gene (location 33,620,000)
		
		TrackFactory.addThickSeparatorTrack(plot.getDataView());
		
		TrackFactory.addPeakTrack(plot, new DataSource(MACS_DATA_FILE));

		TrackFactory.addThickSeparatorTrack(plot.getDataView());

		TrackFactory.addReadTracks(
				plot, 
				Arrays.asList(new DataSource[] { new DataSource(ELAND_DATA_FILE) }),
				Arrays.asList(new DataSource[] { new DataSource(ELAND_DATA_FILE) }),
				new DataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_seq.tsv")
		);
		
		TrackFactory.addRulerTrack(plot);
		plot.start("Y", 1024 * 1024 * 250d);
		plot.moveDataBpRegion(10000L, 10000L);
		
		ChartPanel panel = new ChartPanel(new JFreeChart(plot));
		panel.setPreferredSize(new Dimension(800, 600));
		plot.chartPanel = panel;

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
