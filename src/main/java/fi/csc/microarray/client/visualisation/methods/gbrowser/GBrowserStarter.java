package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedList;

import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JSplitPane;
import javax.swing.WindowConstants;

import org.apache.activemq.kaha.impl.async.DataFile;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Interpretation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Interpretation.TrackType;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.VisualConstants;

/**
 * Quick and dirty starter utility for genome browser development and debugging.
 */
public class GBrowserStarter {

	private static void checkData(DataUrl... urls) throws IOException, URISyntaxException {

		boolean fileNotFoundFail = false;
		for (DataUrl dataUrl : urls) {
			
			File file = dataUrl.getLocalFile();
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
		
		String dataPath = System.getProperty("user.home") + "/example-data/";
		DataUrl BAM_DATA_FILE = new DataUrl(new File(dataPath + "GM12878.bam").toURI().toURL(), "GM12878.bam");
		DataUrl BAI_DATA_FILE = new DataUrl(new File(dataPath + "GM12878.bam.bai").toURI().toURL(), "GM12878.bam.bai");
		DataUrl BED_DATA_FILE = new DataUrl(new File(dataPath + "colors.bed").toURI().toURL(), "colors.bed");
		DataUrl VCF_DATA_FILE = new DataUrl(new File(dataPath + "var.flt.vcf").toURI().toURL(), "var.flt.vcf");
		DataUrl GTF1_DATA_FILE = new DataUrl(new File(dataPath + "cufflinks-gtf/merged-sort.gtf").toURI().toURL(), "merged-sort.gtf");
		//DataUrl GTF2_DATA_FILE = new DataUrl(new File(dataPath + "cufflinks-gtf/transcripts-sort.gtf").toURI().toURL(), "transcripts-sort.gtf");	
		DataUrl GTF2_DATA_FILE = new DataUrl(new File(dataPath + "Homo_sapiens.GRCh37.69-sort.gtf").toURI().toURL().toURI().toURL(), "Homo_sapiens.GRCh37.69-sort.gtf");
		DataUrl GTF3_DATA_FILE = new DataUrl(new URL("http://chipster-filebroker.csc.fi:7060/public/annotations/tmp/Homo_sapiens.GRCh37.69-sort.gtf").toURI().toURL(), "Homo_sapiens.GRCh37.66-sort.gtf");
		 
		DataUrl CNA_DATA_FILE = new DataUrl(new File(dataPath + "cna/regions.tsv").toURI().toURL(), "regions.tsv");
		DataUrl PHENODATA_FILE = new DataUrl(new File(dataPath + "cna/phenodata.tsv").toURI().toURL(), "phenodata.tsv");

		LinkedList<Interpretation> interpretations = new LinkedList<Interpretation>();

		for (int i = 0; i < 2; i++) {
			Interpretation reads = new Interpretation(TrackType.READS, BAM_DATA_FILE);
			reads.setIndexData(BAI_DATA_FILE);			
			interpretations.add(reads);
		}

		//Bed with or without header
		interpretations.add(new Interpretation(TrackType.REGIONS, BED_DATA_FILE));
//		interpretations.add(new Interpretation(TrackType.REGIONS, new DataFile(BED_DATA_FILE)));
//		interpretations.add(new Interpretation(TrackType.REGIONS_WITH_HEADER, new BasicDataFile(data)));		
//		interpretations.add(new Interpretation(TrackType.VCF, new DataFile(VCF_DATA_FILE)));		
//		interpretations.add(new Interpretation(TrackType.GTF, GTF2_DATA_FILE));
//		interpretations.add(new Interpretation(TrackType.GTF, GTF3_DATA_FILE));
//		Interpretation cna = new Interpretation(TrackType.CNA, CNA_DATA_FILE);
//		interpretations.add(cna);

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
