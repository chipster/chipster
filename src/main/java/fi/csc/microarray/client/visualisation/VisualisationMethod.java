package fi.csc.microarray.client.visualisation;

import java.util.List;

import javax.swing.ImageIcon;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.basic.BasicModule.VisualisationMethods;

/**
 * <p>A data visualisation methods. This class is used for
 * showing these as options in the visualisation choice panel, as well as for
 * sending requests to visualiser components and estimating visualisation time.</p>
 * 
 * <p>Method for collecting estimated values for execution time: 1. enable debug outputs in
 * client.visualisation.VisualisationTaskManager 2. perform a lot of
 * operations with different datasets 3. save output to .csv -file and
 * import to excel 4. calculate average time/bytelength value for every
 * visualisation.</p>
 * 
 * @author Janne Käki, Aleksi Kallio, Petri Klemelä
 * 
 */
public class VisualisationMethod {
	
	public String getName() {
		return name;
	}

	private final ClientApplication application = Session.getSession().getApplication();

	private String name;
	private Class<? extends Visualisation> visualiser;
	private String iconPath;
	private ImageIcon icon;
	private int orderNumber;
	/**
	 * Estimated visualisation duration (in milliseconds) per byte of the data length
	 */
	private double durationEstimationFactor;
	private String helpAddress = null;

	public VisualisationMethod(String name, Class<? extends Visualisation> visualiser, String iconPath, int orderNumber, double durationEstimationFactor) {
		this.name = name;
		this.visualiser = visualiser;
		this.iconPath = iconPath;
		this.orderNumber = orderNumber;
		this.durationEstimationFactor = durationEstimationFactor;
	}

	public VisualisationMethod(String name, Class<? extends Visualisation> visualiser, String iconPath, int orderNumber, double durationEstimationFactor, String helpAddress) {
		this(name, visualiser, iconPath, orderNumber, durationEstimationFactor);
		this.helpAddress = helpAddress;
	}
	
	public String toString() {
		return this.name;
	}

	public Class<? extends Visualisation> getVisualiserClass() {
		return visualiser;
	}

	public Visualisation getVisualiser(VisualisationFrame frame) throws Exception {
		
			Visualisation visualisation = visualiser.getConstructor().newInstance();
			visualisation.initialise(frame);
			return visualisation;
	}

	/**
	 * Empty constructor to make the use of visualisations more convenient for
	 * the purposes without direct relation to specific visualisation frame,
	 * e.g. finding out the applicability for particular dataset.
	 * 
	 */
	public Visualisation getHeadlessVisualiser() {
		try {
			return visualiser.getConstructor().newInstance();
		} catch (Exception e) {
			application.reportException(e);
			return null;
		}
	}

	public int getOrderNumber() {
		return orderNumber;
	}

	public ImageIcon getIcon() {
		if (icon == null) {
			icon = VisualConstants.getIcon(iconPath);
		}
		return icon;
	}

	public static VisualisationMethod getDefault() {
		return VisualisationMethods.DATA_DETAILS;
	}
	
	public boolean isDefault() {
		return this == VisualisationMethods.DATA_DETAILS || this == VisualisationMethods.SESSION_DETAILS;
	}

	/**
	 * @param datas
	 * @return long estimated duration of this visualisation in millisecond for
	 *         the given datasets
	 */
	public long estimateDuration(List<DataBean> datas) {
		if (datas.size() > 0) {
			return (long) (Session.getSession().getDataManager().getContentLength(datas.get(0)) * durationEstimationFactor * datas.size());
		}
		return -1;
	}

	public boolean isApplicableTo(DataBean bean) throws MicroarrayException {

		if (bean != null) {
			return this.getHeadlessVisualiser().canVisualise(bean);
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

	public String getHelpAddress() {
		return helpAddress;
	}
}