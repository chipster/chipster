package fi.csc.microarray.client.visualisation;

import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import javax.swing.ImageIcon;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.VisualConstants;
import fi.csc.microarray.client.visualisation.methods.ArrayLayout;
import fi.csc.microarray.client.visualisation.methods.ClusteredProfiles;
import fi.csc.microarray.client.visualisation.methods.Empty;
import fi.csc.microarray.client.visualisation.methods.ExpressionProfile;
import fi.csc.microarray.client.visualisation.methods.HierarchicalClustering;
import fi.csc.microarray.client.visualisation.methods.Histogram;
import fi.csc.microarray.client.visualisation.methods.HtmlViewer;
import fi.csc.microarray.client.visualisation.methods.ImageViewer;
import fi.csc.microarray.client.visualisation.methods.PhenodataEditor;
import fi.csc.microarray.client.visualisation.methods.SOM;
import fi.csc.microarray.client.visualisation.methods.Scatterplot;
import fi.csc.microarray.client.visualisation.methods.Scatterplot3DPCA;
import fi.csc.microarray.client.visualisation.methods.Spreadsheet;
import fi.csc.microarray.client.visualisation.methods.TextViewer;
import fi.csc.microarray.client.visualisation.methods.VennDiagram;
import fi.csc.microarray.client.visualisation.methods.Volcanoplot;
import fi.csc.microarray.client.visualisation.methods.threed.Scatterplot3D;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * An enumeration for all available data visualisation methods. It is used for
 * showing these as options in the visualisation choice panel, as well as for
 * sending requests to visualiser components.
 * 
 * @author Janne KÃ¯Â¿Â½ki, Aleksi Kallio, Petri KlemelÃ¯Â¿Â½
 * 
 */
public enum VisualisationMethod {

	/*
	 * Method for collecting estimated values: 1. enable debug outputs in
	 * client.visualisation.VisualisationTaskManager 2. perform a lot of
	 * operations with different datasets 3. save output to .csv -file and
	 * import to excel 4. calculate average time/bytelength value for every
	 * visualisation
	 */

	// None is always present and spaces reserve space also for the longer names
	// when they aren't present.
	// This is the easiest found way to keep the double visualisation split
	// steady. Duration estimation 0 means unknown.
	NONE("None                                   ", 
			Empty.class, VisualConstants.EMPTY_MENUICON, 0, 0), 
			SPREADSHEET("Spreadsheet", Spreadsheet.class, VisualConstants.SPREADSHEET_MENUICON, 2, 0.0007), 
			PHENODATA("Phenodata editor", PhenodataEditor.class, VisualConstants.PHENODATA_MENUICON, 3, 0), 
			ARRAY_LAYOUT("Array layout", ArrayLayout.class, VisualConstants.ARRAY_MENUICON, -1, 0.0009), 
			HISTOGRAM("Histogram", Histogram.class, VisualConstants.HISTOGRAM_MENUICON, -1, 0.024), 
			SCATTERPLOT("Scatterplot", Scatterplot.class, VisualConstants.SCATTER_MENUICON, -1, 0.039), 
			SCATTERPLOT3D("3D Scatterplot", Scatterplot3D.class, VisualConstants.SCATTER3D_MENUICON, -1, 0.082), 
			SCATTERPLOT3DPCA("3D Scatterplot for PCA", Scatterplot3DPCA.class, VisualConstants.SCATTER3DPCA_MENUICON, -1, 0.082),
			VOLCANOPLOT("Volcano plot", Volcanoplot.class, VisualConstants.VOLCANO_MENUICON, -1, 0.039), 
			SOM("SOM", SOM.class, VisualConstants.SOM_MENUICON, 3, 0.034), 
			HIERARCHICAL("Hierarchical clustering", HierarchicalClustering.class, VisualConstants.HC_MENUICON, 3, 0.09), 
			EXPRESSION_PROFILE("Expression profile", ExpressionProfile.class, VisualConstants.PROFILE_MENUICON, -1, 0.1), 
			CLUSTERED_PROFILES("Clustered profiles", ClusteredProfiles.class, VisualConstants.PROFILES_MENUICON, -1, 0.087), 
			SHOW_IMAGE("Show image", ImageViewer.class, VisualConstants.IMAGE_MENUICON, 1, 0.015), 
			WEBVIEW("View page", HtmlViewer.class, VisualConstants.HTML_MENUICON, 1, 0.008), 
			VIEW_TEXT("View text", TextViewer.class, VisualConstants.TEXT_MENUICON, 1, 0.023), 
			VENN_DIAGRAM("Venn-diagram", VennDiagram.class, VisualConstants.VENN_MENUICON, 1, 0);

	private static LinkedList<VisualisationMethod> orderedDefaultCandidates;

	private final ClientApplication application = Session.getSession().getApplication();

	static {
		orderedDefaultCandidates = new LinkedList<VisualisationMethod>();
		for (VisualisationMethod value : values()) {
			if (value.getOrderNumber() >= 0) {
				orderedDefaultCandidates.add(value);
			}
		}
		Collections.sort(orderedDefaultCandidates, new Comparator<VisualisationMethod>() {
			public int compare(VisualisationMethod method1, VisualisationMethod method2) {
				return new Integer(method2.getOrderNumber()).compareTo(new Integer(method1.getOrderNumber()));
			}
		});
	}

	private String name;
	private Class<? extends Visualisation> visualiser;
	private ImageIcon icon;
	private int orderNumber;
	// Estimated visualisation duration (in milliseconds) per byte of the data
	// length
	private double durationEstimationFactor;

	private VisualisationMethod(String name, Class<? extends Visualisation> visualiser, ImageIcon icon, int orderNumber, double durationEstimationFactor) {
		this.name = name;
		this.visualiser = visualiser;
		this.icon = icon;
		this.orderNumber = orderNumber;
		this.durationEstimationFactor = durationEstimationFactor;
	}

	public String toString() {
		return this.name;
	}

	public Class<? extends Visualisation> getVisualiserClass() {
		return visualiser;
	}

	public Visualisation getVisualiser(VisualisationFrame frame) {
		try {
			return visualiser.getConstructor(VisualisationFrame.class).newInstance(frame);
		} catch (Exception e) {
			application.reportException(e);
			return null;
		}
	}

	/**
	 * Empty constructor to make the use of visualisations more convenient for
	 * the purposes without direct relation to specific visualisation frame,
	 * e.g. finding out the applicability for particular dataset.
	 * 
	 */
	public Visualisation getHeadlessVisualiser() {
		try {
			return visualiser.getConstructor(VisualisationFrame.class).newInstance((Object) null);
		} catch (Exception e) {
			application.reportException(e);
			return null;
		}
	}

	public int getOrderNumber() {
		return orderNumber;
	}

	public ImageIcon getIcon() {
		return icon;
	}

	public static VisualisationMethod getDefault() {
		return VisualisationMethod.NONE;
	}

	/**
	 * @param datas
	 * @return long estimated duration of this visualisation in millisecond for
	 *         the given datasets
	 */
	public long estimateDuration(List<DataBean> datas) {
		if (datas.size() > 0) {
			return (long) (datas.get(0).getContentLength() * durationEstimationFactor * datas.size());
		}
		return -1;
	}

	Map<DataBean, Boolean> canVisualiseCache = new HashMap<DataBean, Boolean>();

	public boolean isApplicableTo(DataBean bean) throws MicroarrayException {

		// The results of the canVisualise are saved to map to speed up
		// the selection of datasets. The maps are emptied only when the client
		// is closed, but this *shouldn't* be a problem for either memory
		// consumption or performance.
		if (bean != null) {
			// Has to be wrapper to allow null
			Boolean result;
			if ((result = canVisualiseCache.get(bean)) == null) {
				result = this.getHeadlessVisualiser().canVisualise(bean);
				canVisualiseCache.put(bean, result);
			}
			return result;
		}
		return false;
	}

	public boolean isApplicableTo(List<DataBean> beans) throws MicroarrayException {
		if (!this.getHeadlessVisualiser().isForMultipleDatas()) {
			return false;
		} else {
			return this.getHeadlessVisualiser().canVisualise(beans);
		}
	}

	public static Iterable<VisualisationMethod> orderedDefaultCandidates() {
		return orderedDefaultCandidates;
	}
	
	public static VisualisationMethod getDefaultVisualisationFor(DataBean dataBean) throws IOException, MicroarrayException {
		for (VisualisationMethod method : VisualisationMethod.orderedDefaultCandidates()) {
			if (method != VisualisationMethod.NONE && method.isApplicableTo(dataBean)) {
				return method;
			}
		}
		return null;
	}
}