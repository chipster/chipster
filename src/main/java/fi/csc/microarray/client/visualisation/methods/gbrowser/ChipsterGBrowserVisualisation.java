package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JPanel;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.selection.IntegratedEntity;
import fi.csc.microarray.client.selection.PointSelectionEvent;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.DataFile;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.Interpretation;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser.TrackType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * Facade class that hides genome browser internals and exposes an API that is compatible 
 * with Chipster visualization system. 
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * @see GenomePlot
 */
public class ChipsterGBrowserVisualisation extends Visualisation {
	
	public static class BeanDataFile extends DataFile {
		
		private DataBean bean;
		
		public BeanDataFile(DataBean data) {
			this.bean = data;
		}

		public String getName() {
			return bean.getName();
		}
		
		/* (non-Javadoc)
		 * Stream has to be closed after use
		 * 
		 * @see fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowser.DataFile#getInputStream()
		 */
		public InputStream getInputStream() throws IOException {
							
			return bean.getContentByteStream();
		}
		public File getLocalFile() throws IOException {
			return Session.getSession().getDataManager().getLocalFile(bean);
		}
	}
	
	public static class DataBeanInterpretation extends Interpretation {

		public List<BeanDataFile> summaryDatas = new LinkedList<BeanDataFile>();

		public DataBeanInterpretation(TrackType type, BeanDataFile primaryData) {
			super(type, primaryData);
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
		protected void reportException(Exception e) {
			application.reportException(e);
		}
		
		@Override
		protected void showDialog(String title, String message, String details, boolean warning, boolean dialogShowDetails, boolean modal) {
			
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
		}
		
		protected void showVisualisation() {

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
		
		protected void runBlockingTask(String taskName, Runnable runnable) {
			application.runBlockingTask(taskName, runnable);
		}
		
		protected void initialiseUserDatas() throws IOException {
			for (Interpretation interpretation : getInterpretations()) {
				initialiseUserData(interpretation.getPrimaryData());
				initialiseUserData(interpretation.getIndexData());
				for (DataFile summaryData : interpretation.getSummaryDatas()) {
					initialiseUserData(summaryData);
				}
			}
		}

		protected void initialiseUserData(DataFile data) throws IOException {
			if (data != null) {				
				data.getLocalFile();
			}
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
			
			if (data.isContentTypeCompatitible("text/plain")) {
				// ELAND result / export
				interpretations.add(new DataBeanInterpretation(TrackType.READS, new BeanDataFile(data)));

			} else if (data.isContentTypeCompatitible("text/bed")) {
				// BED (ChIP-seq peaks)
				interpretations.add(new DataBeanInterpretation(TrackType.REGIONS, new BeanDataFile(data)));

			} else if (data.isContentTypeCompatitible("text/tab")) {
				// peaks (with header in the file)
				interpretations.add(new DataBeanInterpretation(TrackType.REGIONS_WITH_HEADER, new BeanDataFile(data)));

			} else if ((data.isContentTypeCompatitible("application/bam"))) {
				// BAM file
				interpretations.add(new DataBeanInterpretation(TrackType.READS, new BeanDataFile(data)));

			} else if ((data.isContentTypeCompatitible("text/vcf"))) {
				// Vcf file
				interpretations.add(new DataBeanInterpretation(TrackType.VCF, new BeanDataFile(data)));
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
