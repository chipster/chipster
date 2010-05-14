package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.Cursor;
import java.awt.Dimension;
import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;


public class GenomeBrowserStarter {

	private static final File ELAND_DATA_FILE = new File("/home/akallio/Desktop/STAT1/STAT1_treatment_aggregated_filtered_chr1_sorted.txt");
	private static final URL URL_ROOT;

	static {
		try {
			URL_ROOT = new URL("http://chipster-devel.csc.fi:8050/public/annotations");
		} catch (MalformedURLException e) {
			throw new RuntimeException(e);
		}
	}

	public static void main(String[] args) throws IOException {
		GenomePlot plot = new GenomePlot(true);
		TrackFactory.addCytobandTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.57_karyotype.tsv"));
		TrackFactory.addGeneTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_genes.tsv"));
//		TrackFactory.addMirnaTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_miRNA.tsv"));
		TrackFactory.addTranscriptTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_transcripts.tsv"));
//		TrackFactory.addPeakTracks(plot, new DataSource(FILE_ROOT, "results_ar_dht_I_and_II_vs_rigg_combined_p1e5_m10_bw175_ts30_peaks.bed"));
//		TrackFactory.addWigTrack(plot, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_miRNA.fsf"));
//		TrackFactory.addReadTracks(plot, new DataSource(FILE_ROOT, "treatmentdata_FoxA1_sorted_spaced.dat"), new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_seq.tsv"));
//		TrackFactory.addReadTracks(plot, new DataSource(FILE_ROOT, "treatmentdata_FoxA1_sorted_spaced.dat"), new DataSource(FILE_ROOT, "../genomebrowser_data/annotations/Homo_sapiens.GRCh37.56_seq.fsf"));
//		TrackFactory.addReadTracks(plot, new DataSource(FILE_ROOT, "eland_result_sorted.tsv"), 
//				new DataSource(new File("/home/akallio/chipster-share/genomebrowser_data/tsv_annotations/Homo_sapiens.GRCh37.55_seq.tsv")));
		TrackFactory.addReadTracks(plot, new DataSource(ELAND_DATA_FILE), 
				new DataSource(URL_ROOT, "Homo_sapiens.NCBI36.54_seq.tsv"));
		
		TrackFactory.addRulerTrack(plot);
		plot.start("1");
		
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
