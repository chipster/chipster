package fi.csc.microarray.databeans.biobeans;

import java.io.IOException;
import java.util.Arrays;
import java.util.LinkedList;

import org.apache.log4j.Logger;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;

// TODO at some point BioBean should be obsoleted...
public class BioBean {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(BioBean.class);
	
	private static final String CDNA_SAMPLE_FG = "/column/sample";
	private static final String AFFY_INTENSITY = "/column/MEAN";
	public static final String PHENODATA_SAMPLE_COLUMN = "sample";
	public static final String PHENODATA_NAME_COLUMN = "original_name";
	public static final String PHENODATA_DESCRIPTION_COLUMN = "description";
	private static final String PHENODATA_CHIPTYPE_COLUMN = "chiptype";
	private static final String PHENODATA_GROUP_COLUMN = "group";
	
	public static boolean isEditablePhenodataColumn(String columnName) {
		return !PHENODATA_SAMPLE_COLUMN.equals(columnName) && !PHENODATA_NAME_COLUMN.equals(columnName)  && !PHENODATA_CHIPTYPE_COLUMN.equals(columnName);
	}
	
	public static boolean isGroupPhenodataColumn(String columnName) {
		return PHENODATA_GROUP_COLUMN.equals(columnName);
	}
	
	private DataBean dataBean;
	
	public BioBean(DataBean dataBean) {
		this.dataBean = dataBean;
	}
	
	
	public int getColorCount() {
		if (dataBean.queryFeatures(CDNA_SAMPLE_FG).exists()) {
			return 2;
		} else if (dataBean.queryFeatures(AFFY_INTENSITY).exists()) {
			return 1;
		} else {
			return 0;
		}
	}
	
	public VisualisationMethod getDefaultVisualisation() throws IOException, MicroarrayException {
		for (VisualisationMethod method : VisualisationMethod.orderedDefaultCandidates()) {
			if (method != VisualisationMethod.NONE && method.isApplicableTo(dataBean)) {
				return method;
			}
		}
		return null;
	}
	
	/**
	 * Selects the proper source dataset ie. the dataset that is not a hidden phenodata dataset.
	 */
	public DataBean getProperSource() {
		
		logger.debug(dataBean == null ? "databean is null" : ("getting proper source of " + dataBean.getName() + ", source count is " + dataBean.getLinkTargets(Link.DERIVATION).size()));
		
		if (dataBean == null || dataBean.getLinkTargets(Link.DERIVATION).size() == 0) {
			return null;
			
		} else if (dataBean.getLinkTargets(Link.DERIVATION).size() == 1) {
			return dataBean.getLinkTargets(Link.DERIVATION).iterator().next();
			
		} else {
			LinkedList<DataBean> sourceCollector = new LinkedList<DataBean>();
			for (DataBean source : dataBean.getLinkTargets(Link.DERIVATION)) {
				if (source.queryFeatures("/phenodata").exists()) {
					sourceCollector.add(source);
				}
			}
			if (sourceCollector.size() == 0 || sourceCollector.size() > 1) {
				return null; // no definite source was found
			}
			
			return sourceCollector.getFirst(); // return the one and only proper source
		}
	}
	
	/**
	 * Gets the operation history of this dataset as a cronological list
	 * (implemented as vector) of DataSetHistoryWrapper objects, which
	 * practically are DataSets with a slightly altered toString output.
	 * The list will contain all the datasets on the workflow, starting
	 * from raw data and ending at current dataset.
	 * 
	 * @return The chronologically ascending list of dataset history wrappers.
	 */
	public DataBean[] getSourcePath() {
		LinkedList<DataBean> list = new LinkedList<DataBean>();
		if (getProperSource() != null) {
			list.addAll(Arrays.asList(new BioBean(getProperSource()).getSourcePath()));
		}
		list.add(dataBean);
		return list.toArray(new DataBean[0]);
	}	
}
