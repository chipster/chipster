package fi.csc.microarray.databeans.features.bio;

import java.io.IOException;
import java.io.InputStream;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.BoolFalseFeature;
import fi.csc.microarray.databeans.features.BoolTrueFeature;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProviderBase;

/**
 * Checks if bean has embedded binary content i.e. is a text data that has binary
 * content embedded into it. This is for ensuring that we do not try to process such
 * content as text based even if it looks like that.
 * 
 * @author Aleksi Kallio
 *
 */
public class EmbeddedBinaryProvider extends FeatureProviderBase {

	private static final byte[] CEL4_HEADER = {	0x40, 0x00, 0x00, 0x00, 0x04, 0x00, 0x00, 0x00};
	
	public Feature createFeature(String namePostfix, DataBean bean) {
		
		if (!"".equals(namePostfix)) {
			throw new RuntimeException("unknown name postfix: " + namePostfix);
			
		} else {
			
			// check if starts with known embedded binary header
			InputStream contentByteStream = null;
			try {
				contentByteStream = bean.getContentByteStream();
				for (int i = 0; i < CEL4_HEADER.length; i++) {
					// note: EOS => -1 returned from read() => loop terminates
					if (CEL4_HEADER[i] != contentByteStream.read()) {
						return new BoolFalseFeature(bean, this);
					}
				}
				
			} catch (Exception e) {
				throw new RuntimeException(e);
				
			} finally {
				try {
					contentByteStream.close();
				} catch (IOException e) {
					// ignore
				}
			}
			
			// check succeeded, this is embedded binary
			return new BoolTrueFeature(bean, this);
		}
	}

}
