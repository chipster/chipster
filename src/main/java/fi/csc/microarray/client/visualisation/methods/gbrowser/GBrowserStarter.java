package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.LinkedList;

import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JSplitPane;
import javax.swing.WindowConstants;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.DataFile;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.Interpretation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.TrackType;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.VisualConstants;

/**
 * Quick and dirty starter utility for genome browser development and debugging.
 */
public class GBrowserStarter {

	private static void checkData(File... files) {

		boolean fileNotFoundFail = false;
		for (File file : files) {
			if (!file.exists()) {
				System.err.println("File not found: " + file);
				fileNotFoundFail = true;
			}
		}
		if (fileNotFoundFail) {
			System.exit(1);
		}
	}

	public static void main(String[] args) throws Exception {
		
		//Get rid of Chipster logging errors
		DirectoryLayout.initialiseStandaloneClientLayout();

		//Define data
		
		String dataPath = System.getProperty("user.home") + "/chipster/";
		File BAM_DATA_FILE = new File(dataPath + "hg19_chr20.bam");
		File BAI_DATA_FILE = new File(dataPath + "hg19_chr20.bam.bai");
		File VCF_DATA_FILE = new File(dataPath + "var.flt.vcf");
		
//		String dataPath = System.getProperty("user.home") + "/chipster/ohtu/"; 
//		File BAM_DATA_FILE = new File(dataPath + "SRR064438-chr17-chr20.bam");
//		File BAI_DATA_FILE = new File(dataPath + "SRR064438-chr17-chr20.bam.bai");
//		File BED_DATA_FILE = new File(dataPath + "peaks.bed");

		LinkedList<Interpretation> interpretations = new LinkedList<Interpretation>();

		for (int i = 0; i < 2; i++) {
			Interpretation reads = new Interpretation(TrackType.READS, new DataFile(BAM_DATA_FILE));
			reads.setIndexData(new DataFile(BAI_DATA_FILE));			
			interpretations.add(reads);
		}

		//Bed with or without header
//		interpretations.add(new Interpretation(TrackType.REGIONS, new DataFile(BED_DATA_FILE)));
//		interpretations.add(new Interpretation(TrackType.REGIONS, new DataFile(BED_DATA_FILE)));
//		interpretations.add(new BasicInterpretation(TrackType.REGIONS_WITH_HEADER, new BasicDataFile(data)));
		interpretations.add(new Interpretation(TrackType.VCF, new DataFile(VCF_DATA_FILE)));

		checkData(BAM_DATA_FILE, BAI_DATA_FILE);

		//Create gui
		
		JFrame frame = new JFrame();
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);

		BasicGBrowser browser = new BasicGBrowser();
		browser.initialise();

		//Has to be before getVisualisation()
		JComponent parameterPanel =  browser.getParameterPanel();
		JComponent component = browser.getVisualisation(interpretations);

		JSplitPane split = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, component, parameterPanel);
		split.setDividerLocation(1000);
		frame.add(split);
		frame.setSize(1280, 800);
		frame.setVisible(true);
	}

	/**
	 * Configure annotation and icon locations
	 */
	private static class BasicGBrowser extends GBrowser {

		public ImageIcon getIcon(String path) {

			//Chipster packaging specific implementation
			return new ImageIcon(VisualConstants.class.getResource(path));
		}

		public URL getRemoteAnnotationsUrl() throws Exception {
			return new URL("http://chipster-filebroker.csc.fi:8080/public/annotations");
		}

		public File getLocalAnnotationDir() throws IOException {

			//All files in this folder will be DELETED!
			return new File(System.getProperty("user.home") + "/.chipster/annotations");
		}
	}
}
