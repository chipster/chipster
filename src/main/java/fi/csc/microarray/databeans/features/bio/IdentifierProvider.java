package fi.csc.microarray.databeans.features.bio;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProviderBase;
import fi.csc.microarray.databeans.features.NonexistingFeature;

/**
 * <p>
 * Tool for retrieving phenodata related information. Supports following paths:
 * <ul>
 * <li>/ : a valid feature is returned if this bean is identified as phenodata</li>
 * <li>/is-complete : a valid feature is returned if this phenodata is complete</li>
 * <li>/describe/samplename : if phenodata contains descriptions or real sample names it is
 * returned for samplename, otherwise samplename is returned</li>
 * <li>/linked : linked phenodata bean is located and following path is
 * executed against it</li>
 * </ul>
 * </p>
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
