package fi.csc.microarray.databeans.features.bio;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProviderBase;
import fi.csc.microarray.databeans.features.NonexistingFeature;

/**
 * Retrieves the row (gene/probe) identifier column from the dataset.
 *   
 * @author Aleksi Kallio
 * 
 */
public class IdentifierProvider extends FeatureProviderBase {


	private static final String NAMED_IDENTIFIER_COLUMN = "/column/identifier";
	private static final String ROW_NAME_COLUMN = "/column/ ";

	public Feature createFeature(String namePostfix, DataBean bean) {

		if (bean.queryFeatures(ROW_NAME_COLUMN).exists()) {
			return bean.queryFeatures(ROW_NAME_COLUMN).asFeature();
			
		} else if (bean.queryFeatures(NAMED_IDENTIFIER_COLUMN).exists()) {
			return bean.queryFeatures(NAMED_IDENTIFIER_COLUMN).asFeature();
			
		} else {
			return new NonexistingFeature(bean, this);
		}
	}
}
