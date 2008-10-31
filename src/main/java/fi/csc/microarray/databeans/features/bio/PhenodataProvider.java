package fi.csc.microarray.databeans.features.bio;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.visualisation.methods.PhenodataEditor;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.BoolFalseFeature;
import fi.csc.microarray.databeans.features.BoolTrueFeature;
import fi.csc.microarray.databeans.features.ConstantStringFeature;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProviderBase;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.Table;

/**
 * 
 * @author Aleksi Kallio
 * 
 */
public class PhenodataProvider extends FeatureProviderBase {

	private static final String CONVERSION_NAMEFRAGMENT = "_to_";
	private static final String LINKED_PHENODATA_NAME = "linked";

	public Feature createFeature(String namePostfix, DataBean bean) {

		DataBean phenoBean;
		String namePostPostFix;
		if (namePostfix.startsWith(LINKED_PHENODATA_NAME)) {
			phenoBean = LinkUtils.retrieveInherited(bean, Link.ANNOTATION);
			namePostPostFix = namePostfix.substring(LINKED_PHENODATA_NAME.length());
		} else {
			phenoBean = bean;
			namePostPostFix = namePostfix;
		}

		if (namePostfix.contains(CONVERSION_NAMEFRAGMENT)) {
			return createConversionFeature(namePostPostFix, phenoBean);
		} else {
			return createBooleanFeature(namePostPostFix, phenoBean);
		}
	}

	private Feature createBooleanFeature(String namePostfix, DataBean bean) {

		boolean isPhenodata = false;
		boolean isComplete = false;

		// check that data has everything we need
		if (bean.queryFeatures("/column/sample").exists() && bean.queryFeatures("/column/chiptype").exists()) {

			isPhenodata = true;

			try {
				boolean hasEmptyGroups = false;
				QueryResult groupFeature = bean.queryFeatures("/column/group");
				if (groupFeature.exists()) {
					Iterable<String> groups = bean.queryFeatures("/column/group").asStrings();
					for (String group : groups) {
						if ("".equals(group.trim())) {
							hasEmptyGroups = true;
							break;
						}
					}
					isComplete = !hasEmptyGroups;
				}
			} catch (MicroarrayException e) {
				throw new RuntimeException(e);
			}
		}

		boolean returnValue;
		if ("is-complete".equals(namePostfix)) {
			returnValue = isComplete;
		} else {
			returnValue = isPhenodata;
		}

		return returnValue ? new BoolTrueFeature(bean, this) : new BoolFalseFeature(bean, this);
	}

	private Feature createConversionFeature(String namePostfix, DataBean bean) {

		String sampleName = namePostfix.substring(namePostfix.lastIndexOf('/') + 1);
		String originalName = sampleName; // return the plain sample name if
											// nothing was found

		if (bean != null) {
			Table table;
			try {
				table = bean.queryFeatures("/column/*").asTable();

				if (table != null) {
					if (table.hasColumn(PhenodataEditor.PHENODATA_DESCRIPTION_COLUMN)) {
						while (table.nextRow()) {
							if (sampleName.equals(table.getStringValue(PhenodataEditor.PHENODATA_SAMPLE_COLUMN))) {
								originalName = table.getStringValue(PhenodataEditor.PHENODATA_DESCRIPTION_COLUMN);
								break;
							}
						}
					} else if (table.hasColumn(PhenodataEditor.PHENODATA_NAME_COLUMN)) {
						while (table.nextRow()) {
							if (sampleName.equals(table.getStringValue(PhenodataEditor.PHENODATA_SAMPLE_COLUMN))) {
								originalName = table.getStringValue(PhenodataEditor.PHENODATA_NAME_COLUMN);
								break;
							}
						}
					}
				}

			} catch (MicroarrayException e) {
				// do nothing, value was just not found
			}
		}

		return new ConstantStringFeature(bean, this, originalName);
	}

}
