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

	private static final File FILE_ROOT = new File("/home/akallio/chipster-share/genome_browser");
	private static final URL URL_ROOT;

	static {
		try {
			URL_ROOT = new URL("http://chipster-devel.csc.fi:8050/public/space_separated_annotations");
		} catch (MalformedURLException e) {
			throw new RuntimeException(e);
		}
	}

	public static void main(String[] args) throws IOException {
		GenomePlot plot = new GenomePlot(true);
		TrackFactory.addCytobandTracks(plot, new DataSource(URL_ROOT, "cytoband_hg17_sorted.fsf"));
		TrackFactory.addGeneTracks(plot, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_genes.fsf"));
//		TrackFactory.addMirnaTracks(plot.getDataView(), new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_miRNA.fsf"));
//		TrackFactory.addTranscriptTracks(plot.getDataView(), new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_transcripts.fsf"));
//		TrackFactory.addPeakTracks(plot, new DataSource(FILE_ROOT, "results_ar_dht_I_and_II_vs_rigg_combined_p1e5_m10_bw175_ts30_peaks.bed"));
//		TrackFactory.addWigTrack(plot, new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_miRNA.fsf"));
//		TrackFactory.addReadTracks(plot, new DataSource(FILE_ROOT, "treatmentdata_FoxA1_sorted_spaced.dat"), new DataSource(URL_ROOT, "Homo_sapiens.GRCh37.56_seq.fsf"));
		TrackFactory.addReadTracks(plot, new DataSource(FILE_ROOT, "treatmentdata_FoxA1_sorted_spaced.dat"), new DataSource(FILE_ROOT, "../genomebrowser_data/annotations/Homo_sapiens.GRCh37.56_seq.fsf"));
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
