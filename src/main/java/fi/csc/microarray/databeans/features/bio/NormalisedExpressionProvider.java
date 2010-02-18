package fi.csc.microarray.databeans.features.bio;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.CalculatingIterable;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.BasicFeature;
import fi.csc.microarray.databeans.features.FeatureProvider;
import fi.csc.microarray.databeans.features.FeatureProviderBase;
import fi.csc.microarray.databeans.features.QueryResult;
import fi.csc.microarray.databeans.features.CalculatingIterable.CalcOperation;
import fi.csc.microarray.exception.MicroarrayException;

public class NormalisedExpressionProvider extends FeatureProviderBase {
	
	public Feature createFeature(String postfix, DataBean bean) {
		return new NormalisedExpression(bean, this);
	}
	
	public static class NormalisedExpression extends BasicFeature {
		
		private static final String EXPRESSION_COLUMN = "/column/chip.microarray1.cel";

		private static final String RED_CHANNEL_INTENSITY = "/column/sample";
		private static final String RED_CHANNEL_BACKGROUND = "/column/samplebg";
		private static final String GREEN_CHANNEL_INTENSITY = "/column/control";
		private static final String GREEN_CHANNEL_BACKGROUND = "/column/controlbg";
		private static final String AFFY_INTENSITY = "/column/MEAN";

		public NormalisedExpression(DataBean bean, FeatureProvider factory) {
			super(bean, factory);
		}

		@Override
		public Iterable<Float> asFloats() throws MicroarrayException {

			if (getDataBean().queryFeatures(EXPRESSION_COLUMN).exists()) {
				return getDataBean().queryFeatures(EXPRESSION_COLUMN).asFloats();

			} else {
				// non-normalised cDNA, figure out some kind of expression value
				QueryResult rciColumn = getDataBean().queryFeatures(RED_CHANNEL_INTENSITY);
				QueryResult rcbColumn = getDataBean().queryFeatures(RED_CHANNEL_BACKGROUND) ;
				QueryResult gciColumn = getDataBean().queryFeatures(GREEN_CHANNEL_INTENSITY);
				QueryResult gcbColumn = getDataBean().queryFeatures(GREEN_CHANNEL_BACKGROUND) ;
				
				if (rciColumn.exists() && rcbColumn.exists() &&	gciColumn.exists() && gcbColumn.exists()) {
					
					CalculatingIterable redIntensity = new CalculatingIterable( 
							rciColumn.asFloats(), 
							rcbColumn.asFloats(),
							CalcOperation.SUBTRACT); 
					CalculatingIterable  greenIntensity = new CalculatingIterable(
							gciColumn.asFloats(), 
							gcbColumn.asFloats(),
							CalcOperation.SUBTRACT); 
					return new CalculatingIterable(redIntensity, greenIntensity, CalcOperation.SUBTRACT); 

				} else if (getDataBean().queryFeatures(AFFY_INTENSITY).exists()) {
					// non-normalised Affy, return intensity value
					return getDataBean().queryFeatures(AFFY_INTENSITY).asFloats();
					
				} else {
					return null;
				}
			}
		}

	}
}
