package fi.csc.microarray.databeans.features.table;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.ConstantFloatFeature;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProviderBase;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;


/**
 * @author Aleksi Kallio
 */
public class RowCountProvider extends FeatureProviderBase {

	private static final String AT_LEAST_ROWS_CACHENAME = "at-least-rows";

	public Feature createFeature(String namePostfix, DataBean bean) {
		
		try {
			
		if (namePostfix.contains("/")) {
			String[] commands = namePostfix.split("/");
			if ("max".equals(commands[0])) {

				int gtRows = Integer.parseInt(commands[1]);

				Integer cachedCount = (Integer)bean.getFromContentCache(AT_LEAST_ROWS_CACHENAME);
				int rowCount;
				
				if (cachedCount != null && cachedCount >= gtRows) {
					// cached count is good enough
					rowCount = cachedCount; 

				} else {
					// count rows
					Table rowCounter = bean.queryFeatures("/column/*").asTable();
					rowCount = 0;
					while (rowCounter != null && rowCounter.nextRow() && rowCount < gtRows) {
						rowCount++;
					}
					bean.putToContentCache(AT_LEAST_ROWS_CACHENAME, (Integer)rowCount);
				}

				return new ConstantFloatFeature(bean, this, rowCount);
				
			} else {
				throw new IllegalArgumentException("row counting command not understood: " + namePostfix);
			}
		} else {
			throw new IllegalArgumentException("must specify row counting command");
		}

		} catch (MicroarrayException e) {
			throw new RuntimeException(e);
		}
	}
}
