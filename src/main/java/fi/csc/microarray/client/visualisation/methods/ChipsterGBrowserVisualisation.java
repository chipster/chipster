package fi.csc.microarray.client.visualisation.methods;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import javax.swing.ImageIcon;
import javax.swing.JComponent;
import javax.swing.JPanel;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.PluginButton;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.selection.PointSelectionEvent;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.Interpretation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.TrackType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient;

/**
 * Facade class that hides genome browser internals and exposes an API that is compatible 
 * with Chipster visualization system. See class GBrowserStarter for how to start genome browser 
 * outside Chipster. 
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * @see GBrowserPlot
 * @see GBrowserStarter
 */
public class ChipsterGBrowserVisualisation extends Visualisation {
	
	private static final String ANNOTATIONS_PATH = "annotations";
	
	public static class BeanDataFile extends DataUrl {
		
		private DataBean bean;
		
		public BeanDataFile(DataBean data) {
			super(null, data.getName());
			this.bean = data;
		}

		@Override
		public String getName() {
			return bean.getName();
		}
		
		/* (non-Javadoc)
		 * Stream has to be closed after use
		 * 
		 * @see fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowser.DataFile#getInputStream()
		 */
		@Override
		public InputStream getInputStream() throws IOException {
							
			return bean.getContentByteStream();
		}
		
		@Override
		public File getLocalFile() throws IOException {
			return Session.getSession().getDataManager().getLocalFile(bean);
		}
		
		@Override
		public URL getUrl() {
			try {
				return getLocalFile().toURI().toURL();
			} catch (MalformedURLException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return null;
		}
	}
	
	public static class DataBeanInterpretation extends Interpretation {

		public List<BeanDataFile> summaryDatas = new LinkedList<BeanDataFile>();

		public DataBeanInterpretation(TrackType type, BeanDataFile beanDataFile) {
			super(type, beanDataFile);
		}
	}
	
	/**
	 * GenomeBrowser should use Chipster's general components when its run in Chipster, but shouldn't 
	 * depend on those when it's run alone. The general version of this functionality is implemented in 
	 * class GenomeBrowser and those implementations are overridden here with Chipster specific calls.   
	 * 
	 * @author klemela
	 */
	private static class ChipsterGBrowser extends GBrowser implements PropertyChangeListener {
		
		private ClientApplication application;

		public ChipsterGBrowser() {
			this.application = Session.getSession().getApplication();
		}
		
		@Override
		public void reportException(Exception e) {
			application.reportException(e);
		}
		
		@Override
		public void showDialog(String title, String message, String details, boolean warning, boolean dialogShowDetails, boolean modal, boolean closeBrowser) {
			
			Severity severity;
			
			if (warning) {
				severity = Severity.WARNING;
			} else {
				severity = Severity.INFO;
			}
			
			DetailsVisibility detailsVisibility;
			
			if (dialogShowDetails) {
				detailsVisibility = DetailsVisibility.DETAILS_VISIBLE;
			} else {
				detailsVisibility = DetailsVisibility.DETAILS_ALWAYS_HIDDEN;
			}
			
			application.showDialog(title, message, details, severity, modal, detailsVisibility, null);
			
			if (closeBrowser) {
				application.setVisualisationMethod(VisualisationMethod.NONE, null, application.getSelectionManager().getSelectedDataBeans(), FrameType.MAIN);
			}
		}
		
		public void showVisualisation() {

			super.showVisualisation();
							
			// Add selection listener (but try to remove first old one that would prevent garage collection of the visualization) 
			application.removeClientEventListener(this);
			application.addClientEventListener(this);
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
					
					setLocation(chr, start, end);
				}
			}
		}
		
		public void runBlockingTask(String taskName, Runnable runnable) {
			application.runBlockingTask(taskName, runnable);
		}
		
		public void initialiseUserDatas() throws IOException {
			for (Interpretation interpretation : getInterpretations()) {
				initialiseUserData(interpretation.getPrimaryData());
				initialiseUserData(interpretation.getIndexData());
				for (DataUrl summaryData : interpretation.getSummaryDatas()) {
					initialiseUserData(summaryData);
				}
			}
		}

		protected void initialiseUserData(DataUrl data) throws IOException {
			if (data != null) {				
				try {
					data.getLocalFile();
				} catch (URISyntaxException e) {
					e.printStackTrace();
				}
			}
		}
		
		public ImageIcon getIcon(String path) {
			return new ImageIcon(VisualConstants.class.getResource(path));
		}
		
		public void openDownloadAnnotationsDialog(final Genome genome) {
			Session.getSession().getApplication().showDialog(
					"Download annotations for " + genome + "?",
					"Downloading annotations is highly recommended to get optimal performace with genome browser.\n\nYou only need to download annotations once, after that they are stored on your local computer for further use.",
					"", Severity.INFO, true, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, new PluginButton() {

						@Override
						public void actionPerformed() {
							try {
								getAnnotationManager().downloadAnnotations(genome);
							} catch (IOException e) {
								throw new RuntimeException(e);
							}
						}

						@Override
						public String getText() {
							return "Download ";
						}
					});

		}
		
		public URL getRemoteAnnotationsUrl() throws Exception {
			FileBrokerClient fileBroker = Session.getSession().getServiceAccessor().getFileBrokerClient();
			
			URL publicURL = fileBroker.getPublicUrl();
			if (publicURL != null) {
				return new URL(publicURL + "/" + ANNOTATIONS_PATH);

			} else {
				return null;
			}
		}
		
		public File getLocalAnnotationDir() throws IOException {
			return DirectoryLayout.getInstance().getLocalAnnotationDir();
		}
	}
	
	private ChipsterGBrowser browser;

	private final ClientApplication application = Session.getSession()
			.getApplication();


	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);

		browser = new ChipsterGBrowser();
		browser.initialise();
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		return getVisualisation(Arrays.asList(new DataBean[] { data }));
	}


	@Override
	public JComponent getVisualisation(java.util.List<DataBean> datas) throws Exception {
		
		return browser.getVisualisation(interpretUserDatas(datas));
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
			
			if (data.isContentTypeCompatitible("text/bed")) {
				// BED (ChIP-seq peaks)
				interpretations.add(new DataBeanInterpretation(TrackType.REGIONS, new BeanDataFile(data)));

			} else if ((data.isContentTypeCompatitible("application/bam"))) {
				// BAM file
				interpretations.add(new DataBeanInterpretation(TrackType.READS, new BeanDataFile(data)));

			} else if ((data.isContentTypeCompatitible("text/vcf"))) {
				// Vcf file
				interpretations.add(new DataBeanInterpretation(TrackType.VCF, new BeanDataFile(data)));
			} else if ((data.isContentTypeCompatitible("text/gtf"))) {
				// Gtf file
				interpretations.add(new DataBeanInterpretation(TrackType.GTF, new BeanDataFile(data)));
			}						
		}

		// Find interpretations for all secondary data types
		for (DataBean data : datas) {

			// Find the interpretation to add this secondary data to
			Interpretation primaryInterpretation = null;
			for (Interpretation interpretation : interpretations) {
				if (data.getName().startsWith(interpretation.getPrimaryData().getName())) {
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
				primaryInterpretation.getSummaryDatas().add(new BeanDataFile(data));

			} else if ((data.isContentTypeCompatitible("application/octet-stream")) &&
					(isIndexData(data))) {
				// BAI file
				if (primaryInterpretation.getIndexData() != null) {
					return null; // already taken, could not bind this secondary data to any primary data
				}
				primaryInterpretation.setIndexData(new BeanDataFile(data));
			}
		}

		// Check that interpretations are now fully initialised
		for (Interpretation interpretation : interpretations) {
			if (interpretation.getPrimaryData().getName().endsWith(".bam") && interpretation.getIndexData() == null) {
				
				String indexName = interpretation.getPrimaryData().getName().replace(".bam", ".bam.bai");
				DataBean indexBean = application.getDataManager().getDataBean(indexName);
				
				if (indexBean == null) {
				
					return null; // BAM is missing BAI
				} else {
					interpretation.setIndexData(new BeanDataFile(indexBean));
				}
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
	public void removeVisualisation() {

		super.removeVisualisation();

		browser.removeVisualisation();

		application.removeClientEventListener(browser);

	}
	
	@Override
	public JPanel getParameterPanel() {
		return browser.getParameterPanel();
	}
}
