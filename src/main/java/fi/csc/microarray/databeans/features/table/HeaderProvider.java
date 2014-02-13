package fi.csc.microarray.databeans.features.table;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.features.ConstantStringFeature;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProviderBase;
import fi.csc.microarray.databeans.features.table.TableColumnProvider.MatrixParseSettings;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.LookaheadLineReader;


public class HeaderProvider extends FeatureProviderBase {
	/**
	 * Logger for this class
	 */
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(HeaderProvider.class);

	public Feature createFeature(String namePostfix, DataBean bean) {
			
		try (BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(Session.getSession().getDataManager().getContentStream(bean, DataNotAvailableHandling.EMPTY_ON_NA)))) {
			
			MatrixParseSettings settings = TableColumnProvider.inferSettings(bean);
			LookaheadLineReader source = new LookaheadLineReader(bufferedReader);
			String header = TableColumnProvider.getHeader(source, settings);
			
			return new ConstantStringFeature(bean, this, header);
			
		} catch (MicroarrayException | IOException e) {
			throw new RuntimeException(e.getMessage() + " (when reading data header out of " + bean.getName() + ")", e);
		}		
	}
}
