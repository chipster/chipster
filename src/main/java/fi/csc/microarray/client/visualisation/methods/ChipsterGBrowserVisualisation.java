package fi.csc.microarray.client.visualisation.methods;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Arrays;
import java.util.HashSet;
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
import fi.csc.microarray.client.selection.DataSelectionManager;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.selection.IntegratedSelectionManager;
import fi.csc.microarray.client.selection.PointSelectionEvent;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowserStarter;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.AnnotationManager.Genome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.BrowserSelectionListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Interpretation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Interpretation.TrackType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Selectable;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.module.chipster.MicroarrayModule;

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
		
		public BeanDataFile(DataBean data, String name) {
			super(null, name);
			this.bean = data;
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
			//Chipster2 backport fix
			return Session.getSession().getDataManager().getLocalFile(bean);
		}
		
		@Override
		public URL getUrl() throws IOException {
			//Chipster2 backport fix
			return getLocalFile().toURI().toURL();
		}

		public DataBean getDataBean() {
			return bean;
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
	private static class ChipsterGBrowser extends GBrowser implements PropertyChangeListener, BrowserSelectionListener {
		
		private ClientApplication application;
		private List<DataBean> datas;

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
							
			// Add selection listener (but try to remove first old one that would prevent garbage collection of the visualization) 
			application.removeClientEventListener(this);
			application.addClientEventListener(this);
			
			getSelectionManager().addSelectionListener(this);
			
			//Update selections for every data
			for (DataBean bean : datas) {
				updateSelectionsFromChipster(bean, null);
			}
		}
		
		@Override
		public void selectionChanged(DataUrl data, Selectable changedSelectable, Object eventSource) {
		
			//Selection changed in genome browser and notify Chipster about it
			
			if (eventSource != this) {
				//The event came from Chipster. Do not send it forward.
				return;
			}
			
			DataSelectionManager dataSelectionManager = Session.getSession().getApplication().getSelectionManager();
						
			DataBean dataBean = null;

			if (data instanceof BeanDataFile) {
				BeanDataFile dataBeanFile = (BeanDataFile) data;
				dataBean = dataBeanFile.getDataBean();
			} else {
				return;
			}

			IntegratedSelectionManager rowSelectionManager = dataSelectionManager.getSelectionManager(dataBean);

			HashSet<Integer> selected = new HashSet<>();								
			for (Selectable selectable : getSelectionManager().getSelectableSet(data)) {

				Integer row = selectable.getIndexKey().getRowNumber();

				//Row number is available only with small files (using InMemoryIndex)
				if (row != null) {					
					selected.add(row);
				}
			}			
			rowSelectionManager.setSelected(selected, eventSource);			
		}			
		
		@Override
		public void propertyChange(PropertyChangeEvent event) {
			
			//Do not process this event, if it originated from the genome browser
			if (event instanceof PointSelectionEvent && event.getSource() != this) {

				PointSelectionEvent pse = (PointSelectionEvent) event;
				DataBean bean = pse.getData();
				
				updateSelectionsFromChipster(bean, event.getSource());
			}		
		}

		private void updateSelectionsFromChipster(DataBean bean, Object eventSource) {
			
			/*
			 * Update selections
			 */
			
			DataSelectionManager dataSelectionManager = application.getSelectionManager();
			IntegratedSelectionManager rowSelectionManager = dataSelectionManager.getSelectionManager(bean);
			
			BeanDataFile data = new BeanDataFile(bean);
			
			//Convert int array to Integer set
			HashSet<Integer> rows = new HashSet<Integer>();
			for (int i : rowSelectionManager.getSelectionAsRows()) {					
				rows.add((Integer)i);
			}								

			getSelectionManager().setRowSelections(data, rows, eventSource);
							
			/*
			 * Not a selection but a request to move
			 */
			
			IntegratedEntity sel = application.getSelectionManager().getSelectionManager(null).getPointSelection();

			if (sel != null && sel.containsKey("chromosome") && sel.containsKey("start")) {

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
		
		public void runBlockingTask(String taskName, Runnable runnable) {
			application.runBlockingTask(taskName, runnable);
		}
		
		public void initialiseUserDatas() throws IOException {
			for (Interpretation interpretation : getInterpretations()) {
				initialiseUserData(interpretation.getPrimaryData());
				initialiseUserData(interpretation.getIndexData());
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
		
		@Deprecated
		public URL getRemoteAnnotationsUrl() throws Exception {
			FileBrokerClient fileBroker = Session.getSession().getServiceAccessor().getFileBrokerClient();
			
			List<URL> publicFiles = fileBroker.getPublicFiles();
			if (publicFiles != null) {
				
				//find only the annotations folder for now
				for (URL url : publicFiles) {
					if  (url.getPath().contains("/" + ANNOTATIONS_PATH)) {
						
						String urlString = url.toString();
						String annotationString = urlString.substring(0, urlString.indexOf("/" + ANNOTATIONS_PATH) + ANNOTATIONS_PATH.length() + 1);
						return new URL(annotationString);
					}
				}
			}
			
			return null;			
		}

		public List<URL> getRemoteAnnotationFiles() throws Exception {
			FileBrokerClient fileBroker = Session.getSession().getServiceAccessor().getFileBrokerClient();
			
			return fileBroker.getPublicFiles();			
		}
		
		
		public File getLocalAnnotationDir() throws IOException {
			return DirectoryLayout.getInstance().getLocalAnnotationDir();
		}
		
		@Override
		public LinkedList<String> getSampleNames(LinkedList<String> sampleNames, DataUrl dataUrl) {
						
			DataBean bean = ((BeanDataFile)dataUrl).bean;
			
			try {
				for (int i = 0; i < sampleNames.size(); i++) {

					String internalName = sampleNames.get(i);

					String sampleName;
					sampleName = bean.queryFeatures("/phenodata/linked/describe/" + internalName).asString();

					sampleNames.set(i, sampleName);				
				}

			} catch (MicroarrayException e) {
				//Use internal names
			}
			
			return sampleNames;
		}

		public void setDatas(List<DataBean> datas) {
			this.datas = datas;
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

		browser.setDatas(datas);
		
		List<Interpretation> interpretations = interpretUserDatas(datas);
		
		if (interpretations != null) {
			return browser.getVisualisation(interpretations);
		} else {
			return null;
		}
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
											
			} else if (data.isContentTypeCompatitible("text/tab")) {
					
				if (data.hasTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_FIRST_TABLE_COLUMN) &&
					data.hasTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_SECOND_TABLE_COLUMN) &&
					data.hasTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_THIRD_TABLE_COLUMN)) {
						
					interpretations.add(new DataBeanInterpretation(TrackType.TSV, new BeanDataFile(data)));	
				}
							
				if (data.hasTypeTag(MicroarrayModule.TypeTags.CHROMOSOME_IN_SECOND_TABLE_COLUMN) &&
					data.hasTypeTag(MicroarrayModule.TypeTags.START_POSITION_IN_THIRD_TABLE_COLUMN) &&
					data.hasTypeTag(MicroarrayModule.TypeTags.END_POSITION_IN_FOURTH_TABLE_COLUMN)) {
					
					interpretations.add(new DataBeanInterpretation(TrackType.TSV_WITH_ROW_ID, new BeanDataFile(data)));
				}
			
				if (data.hasTypeTag(MicroarrayModule.TypeTags.CNA)) {

					// Cna file				
					interpretations.add(new DataBeanInterpretation(TrackType.CNA, new BeanDataFile(data, data.getName())));
				}
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
				
				//Chipster2 backport fixes on following 20 lines
				LinkedList<DataBean> beanList = application.getDataManager().getDataBeans(indexName);
				
				if (beanList.size() == 1) {
					
					DataBean indexBean = beanList.get(0);
					interpretation.setIndexData(new BeanDataFile(indexBean));
					
				} else if (beanList.size() > 1) { 					
					if (browser != null) {						
						//A real visualization attempt, not just applicability check
						browser.showDialog(
								"Unable to determine index file"  , 
								"There are several index files with name '" + indexName + "'. " +
										"Please identify the right index file by selecting it or rename bam and bai file pairs with unique names." , 
										null, false, false, true, true);
						return null;
					}
					
				} else {
					if (browser != null) {						
						//A real visualization attempt, not just applicability check

						browser.showDialog("Missing index file", 
								"There is no index file for data '" + interpretation.getPrimaryData().getName() + "'.",
								null, false, false, true, true);
					}
					return null;
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
