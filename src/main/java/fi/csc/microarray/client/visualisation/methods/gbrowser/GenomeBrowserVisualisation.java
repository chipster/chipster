package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.awt.CardLayout;
import java.awt.Cursor;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JScrollPane;

import org.jfree.chart.JFreeChart;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.selection.PointSelectionEvent;
import fi.csc.microarray.client.visualisation.NonScalableChartPanel;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowser.Interpretation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowser.TrackType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.AnnotationType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.BrowserLauncher;

/**
 * Facade class that hides genome browser internals and exposes an API that is compatible 
 * with Chipster visualization system. 
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * @see GenomePlot
 */
public class GenomeBrowserVisualisation extends Visualisation implements PropertyChangeListener {
	
	public static class DataBeanInterpretation extends Interpretation {

		public List<DataBean> summaryDatas = new LinkedList<DataBean>();
		public DataBean primaryData;
		public DataBean indexData;

		public DataBeanInterpretation(TrackType type, DataBean primaryData) {
			super(type, null);
			this.primaryData = primaryData;
		}

	}
	
	private class ChipsterGenomeBrowser extends GenomeBrowser {
		@Override
		protected void reportException(Exception e) {
			application.reportException(e);
		}
		
		@Override
		protected void showDialog(String title, String message, String details, boolean showDetails, boolean modal) {
			
			if (showDetails) {
				application.showDialog(title, message, details, Severity.WARNING, modal);
			} else {
				application.showDialog(title, message, details, Severity.WARNING, modal, DetailsVisibility.DETAILS_HIDDEN, null);
			}
		}
		
		@Override
		protected void openExternalBrowser(String url) {

			try {
				BrowserLauncher.openURL(url);
			} catch (Exception e) {
				application.reportException(e);
			}
		}
		
		protected void showVisualisation() {

			super.showVisualisation();
							
			// Add selection listener (but try to remove first old one that would prevent garage collection of the visualization) 
			application.removeClientEventListener(GenomeBrowserVisualisation.this);
			application.addClientEventListener(GenomeBrowserVisualisation.this);
		}
		
		protected void runBlockingTask(String taskName, Runnable runnable) {
			application.runBlockingTask(taskName, runnable);
		}

	}
	
	private ChipsterGenomeBrowser browser;

	private final ClientApplication application = Session.getSession()
			.getApplication();


	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);

		browser.initialise();
	}

	
	public void openExternalBrowser(AnnotationType browser) {

		String url = getExternalLinkUrl(browser);	
		Region region = this.plot.getDataView().getBpRegion();
		url = url.replace(AnnotationManager.CHR_LOCATION, region.start.chr.toNormalisedString());
		url = url.replace(AnnotationManager.START_LOCATION, region.start.bp.toString());
		url = url.replace(AnnotationManager.END_LOCATION, region.end.bp.toString());

		try {
			BrowserLauncher.openURL(url);
		} catch (Exception e) {
			application.reportException(e);
		}
	}


	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		return getVisualisation(Arrays.asList(new DataBean[] { data }));
	}


	@Override
	public JComponent getVisualisation(java.util.List<DataBean> datas) throws Exception {
		
		return browser.getVisualisation(interpretUserDatas(datas));
	}


	private void initialiseUserDatas() throws IOException {
		for (Interpretation interpretation : interpretations) {
			initialiseUserData(interpretation.primaryData);
			initialiseUserData(interpretation.indexData);
			for (DataBean summaryData : interpretation.summaryDatas) {
				initialiseUserData(summaryData);
			}
		}
	}

	private void initialiseUserData(DataBean data) throws IOException {
		if (data != null) {
			Session.getSession().getDataManager().getLocalFile(data);
		}
	}


	/**
	 * Create DataSource either for SAM/BAM or ELAND data files.
	 * 
	 * @param tracks
	 * 
	 * @param file
	 * @return
	 * @throws MicroarrayException
	 *             if index file is not selected properly
	 * @throws IOException
	 *             if opening data files fails
	 * @throws URISyntaxException 
	 */
	public DataSource createReadDataSource(DataBean data, DataBean indexData, List<Track> tracks)
			throws MicroarrayException, IOException, URISyntaxException {
		DataSource dataSource = null;

		// Convert data bean into file
		File file = data == null ? null : Session.getSession().getDataManager().getLocalFile(data);

		URL fileUrl = file.toURI().toURL();

		if (data.getName().contains(".bam-summary")) {
			dataSource = new TabixSummaryDataSource(fileUrl);

		} else if (data.getName().contains(".bam") || data.getName().contains(".sam")) {
			File indexFile = Session.getSession().getDataManager().getLocalFile(indexData);
			URL indexFileUrl = indexFile.toURI().toURL();
			dataSource = new SAMDataSource(fileUrl, indexFileUrl);

		} else {
			dataSource = new ChunkDataSource(fileUrl, new ElandParser(), ChunkTreeHandlerThread.class);
		}

		return dataSource;
	}

	private boolean isIndexData(DataBean bean) {
		return bean.getName().endsWith(".bai");
	}

	@Override
	public boolean canVisualise(DataBean data) throws MicroarrayException {
		return canVisualise(Arrays.asList(new DataBean[] { data }));
	}

	@Override
	public boolean canVisualise(java.util.List<DataBean> datas) throws MicroarrayException {
		return interpretUserDatas(datas) != null;
	}

	public class ObjVariable extends Variable {

		public Object obj;

		public ObjVariable(Object obj) {
			super(null, null);
			this.obj = obj;
		}
	}

	private List<Interpretation> interpretUserDatas(List<DataBean> datas) {
		LinkedList<Interpretation> interpretations = new LinkedList<Interpretation>();

		// Find interpretations for all primary data types
		for (DataBean data : datas) {
			


			if (data.isContentTypeCompatitible("text/plain")) {
				// ELAND result / export
				interpretations.add(new DataBeanInterpretation(TrackType.READS, data));

			} else if (data.isContentTypeCompatitible("text/bed")) {
				// BED (ChIP-seq peaks)
				interpretations.add(new DataBeanInterpretation(TrackType.REGIONS, data));

			} else if (data.isContentTypeCompatitible("text/tab")) {
				// peaks (with header in the file)
				interpretations.add(new DataBeanInterpretation(TrackType.REGIONS_WITH_HEADER, data));

			} else if ((data.isContentTypeCompatitible("application/bam"))) {
				// BAM file
				interpretations.add(new DataBeanInterpretation(TrackType.READS, data));

			} else if ((data.isContentTypeCompatitible("text/vcf"))) {
				// Vcf file
				interpretations.add(new DataBeanInterpretation(TrackType.VCF, data));
			}
		}

		// Find interpretations for all secondary data types
		for (DataBean data : datas) {

			// Find the interpretation to add this secondary data to
			Interpretation primaryInterpretation = null;
			for (Interpretation interpretation : interpretations) {
				if (data.getName().startsWith(interpretation.primaryData.getName())) {
					primaryInterpretation = interpretation;
					break;
				}
			}

			if (primaryInterpretation == null) {
				return null; // could not bound this secondary data to any primary data
			}

			if ((data.isContentTypeCompatitible("application/octet-stream")) &&
					(data.getName().contains(".bam-summary"))) {
				// BAM summary file (from custom preprocessor)
				primaryInterpretation.summaryDatas.add(data);

			} else if ((data.isContentTypeCompatitible("application/octet-stream")) &&
					(isIndexData(data))) {
				// BAI file
				if (primaryInterpretation.indexData != null) {
					return null; // already taken, could not bind this secondary data to any primary data
				}
				primaryInterpretation.indexData = data;
			}
		}

		// Check that interpretations are now fully initialised
		for (Interpretation interpretation : interpretations) {
			if (interpretation.primaryData.getName().endsWith(".bam") && interpretation.indexData == null) {
				return null; // BAM is missing BAI
			}
		}


		return interpretations;
	}

	@Override
	public boolean isForSingleData() {
		return true;
	}

	@Override
	public boolean isForMultipleDatas() {
		return true;
	}

	@Override
	public void propertyChange(PropertyChangeEvent event) {
		if (event instanceof PointSelectionEvent) {

			IntegratedEntity sel = application.getSelectionManager().getSelectionManager(null).getPointSelection();

			// Check if we can process this
			if (sel.containsKey("chromosome") && sel.containsKey("start")) {

				Chromosome chr = new Chromosome(sel.get("chromosome"));
				Long start = Long.parseLong(sel.get("start"));

				Long end;
				if (sel.containsKey("end")) {
					end = Long.parseLong(sel.get("end"));
				} else {
					end = null;
				}
				
				browser.setLocation(chr, start, end);
			}
		}
	}

	@Override
	public void removeVisualisation() {

		super.removeVisualisation();

		browser.removeVisualisation();

		application.removeClientEventListener(this);

	}
}
