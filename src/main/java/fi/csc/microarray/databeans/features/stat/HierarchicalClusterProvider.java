package fi.csc.microarray.databeans.features.stat;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.BasicFeature;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProviderBase;
import fi.csc.microarray.databeans.features.NonexistingFeature;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.Strings;

public class HierarchicalClusterProvider extends FeatureProviderBase {

	public Feature createFeature(String namePostfix, final DataBean bean) {
		
		if ("heatmap".equals(namePostfix)) {
			DataBean source = MicroarrayModule.getProperSource(bean);
			if (source != null) {
				return source.queryFeatures("/column/*").asFeature();
			} else {
				return new NonexistingFeature(bean, this);
			}
			
		} else if ("tree".equals(namePostfix)) {
			return new BasicFeature(bean, this) {
				public Iterable<String> asStrings() throws MicroarrayException {
					
					BufferedReader reader = null;
					String tree = "";
					try {
						reader = new BufferedReader(new InputStreamReader(bean.getContentByteStream()));
						boolean first = true;
						for (String line = reader.readLine(); line != null; line = reader.readLine()) {
							tree += line;
							if (first) {
								// check first line to verify that this is a tree - a more robust yet efficient check would be very difficult to do...
								if (!line.contains(":") && !line.contains("(")) {
									throw new MicroarrayException("cannot parse tree: " + Strings.crop(line, 20));
								}
							}
							first = false;
						}						
					} catch (IOException e) {
						throw new MicroarrayException(e);
					} finally {
						try {
							reader.close();
						} catch (Exception e) {
							// ignore
						}
					}
					return Arrays.asList(new String[] {tree});
				}
			};

		} else {
			throw new IllegalArgumentException("unknown postfix: " + namePostfix); 
		}
	}

}
