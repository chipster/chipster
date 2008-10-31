package fi.csc.microarray.databeans.features.stat;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.ConstantTableFeature;
import fi.csc.microarray.databeans.features.Feature;
import fi.csc.microarray.databeans.features.FeatureProviderBase;
import fi.csc.microarray.databeans.features.NonexistingFeature;

public class SomClusterProvider extends FeatureProviderBase {

	private static class Cluster {
		String colour;
		int number;
		float distanceToFirst;
		String values = "";
	}
	
	public Feature createFeature(String namePostfix, DataBean bean) {
		// check that data has everything we need
		if (!bean.queryFeatures("/column/colours").exists() || !bean.queryFeatures("/column/distance2first").exists() ||
				!bean.queryFeatures("/column/cluster").exists() || !bean.queryFeatures("/identifier").exists() ||
				!bean.queryFeatures("/column/griddim").exists()) {
			return new NonexistingFeature(bean, this); // SOM feature is not supported
		}
		
		try {
			
			HashMap<Integer, Cluster> clusterMap = new HashMap<Integer, Cluster>(); 

			// read grid dimensions
			Iterator<Float> dimensions = bean.queryFeatures("/column/griddim").asFloats().iterator();
			int width = dimensions.next().intValue();
			int height = dimensions.next().intValue();
			
			// read clusters
			Iterator<String> colours = bean.queryFeatures("/column/colours").asStrings().iterator();
			Iterator<Float> distances = bean.queryFeatures("/column/distance2first").asFloats().iterator();
			
			for (int i = 1; ; i++) {
				String colour = "";
				if (colours.hasNext()) {
					colour = colours.next();
				}
				if ("".equals(colour.trim()) || "NaN".equals(colour)) {
					// FIXME we should not be seeing NaN's here (but we are...) 
					break; // we break at end-of-iteration or an empty string
				}

				Cluster cluster = new Cluster();
				cluster.number = i;
				cluster.colour = colour;
				cluster.distanceToFirst = distances.next();
				clusterMap.put(i, cluster);
			}
			
			// place genes into clusters
			Iterable<Float> clusters = bean.queryFeatures("/column/cluster").asFloats();
			Iterator<String> genes = bean.queryFeatures("/identifier").asStrings().iterator();

			for (Float clusterNumber : clusters) {
				int number = clusterNumber.intValue();
				Cluster cluster = clusterMap.get(number);
				if (cluster == null) {
					throw new RuntimeException("illegal SOM dataset \"" + bean.getName() + "\": cluster " + number + " was referenced, but does not exist");
				}
				cluster.values += (" " + genes.next()); 
			}
			
			// create data
			LinkedList<Integer> xColumn = new LinkedList<Integer>();
			LinkedList<Integer> yColumn = new LinkedList<Integer>();
			LinkedList<String> colorColumn = new LinkedList<String>();
			LinkedList<String> vectorColumn = new LinkedList<String>();
			LinkedList<String> valuesColumn = new LinkedList<String>();
			
			Iterator<Cluster> clusterIter = clusterMap.values().iterator();
			for (int x = 1; x < width+1; x++) {
				for (int y = 1; y < height+1; y++) {
					Cluster cluster = clusterIter.next();
					xColumn.add(x);
					yColumn.add(y);
					colorColumn.add(cluster.colour);
					vectorColumn.add(""+cluster.distanceToFirst);
					valuesColumn.add(cluster.values);
				}
			}			
			
			// write data
			String[] columns = new String[] { "x", "y", "color", "vector", "values"};
			Object[][] values = new Object[][] { xColumn.toArray(new Integer[0]),
					yColumn.toArray(new Integer[0]), colorColumn.toArray(new String[0]),
					vectorColumn.toArray(new String[0]), valuesColumn.toArray(new String[0]) }; 
			
			return new ConstantTableFeature(bean, this, columns, values);

		} catch (MicroarrayException e) {
			throw new RuntimeException(e);
		}
	}

}
