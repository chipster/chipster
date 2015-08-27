package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.util.ArrayList;

import javax.swing.JComponent;

import org.jfree.chart.BioChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.SOMPlot;
import org.jfree.data.som.SOMDataItem;
import org.jfree.data.som.SOMDataset;

import fi.csc.microarray.client.visualisation.TableAnnotationProvider;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class SOM extends Visualisation {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	/**
	 * Makes a chartPanel with SOM visualization from the data of dataBean.
	 * Finds the biggest coordinate pair from the imput and assumes that all lower
	 * coordinate pairs are found from the data ( otherwise Viski-library throws nullPointer
	 * from the deSelectAll-method of the SOMDataSet class).
	 */
	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {

		// iterate over input data to find out the dimensions of the table
		int maxX=0;
		int maxY=0;
		try (Table som = data.queryFeatures("/clusters/som").asTable()) {

			while (som.nextRow()) {
				if (som.getIntValue("x") > maxX) {
					maxX = som.getIntValue("x");
				}
				if (som.getIntValue("y") > maxY) {
					maxY = som.getIntValue("y");
				}
			}
		}

		// these has to be exactly right, otherwise NPE is thrown
		SOMDataset dataset = new SOMDataset(maxX,maxY);

		// iterate again to read actual data
		try (Table som = data.queryFeatures("/clusters/som").asTable()) {			

			TableAnnotationProvider annotationProvider = new TableAnnotationProvider(data);

			// parse and convert values from dataBean to SOMDataset
			while (som.nextRow()) {

				// color from hex to Color-object
				String colorStr = som.getStringValue("color");			
				Color color;

				if (colorStr.charAt(0) == '#'){
					color = Color.decode(colorStr.substring(0));
				} else {
					throw new RuntimeException("color format not supported for SOM visualization: " + colorStr);
				}			

				// values and vectors from string to typed tables

				ArrayList<String> values = new ArrayList<String>();
				//For all clusters
				for(String idList: som.getStringValue("values").split(", ")){				
					//For all genes in a cluster
					for(String id: idList.trim().split(" ")){
						values.add(annotationProvider.getAnnotatedRowname(id.trim()));
					}
				}
				String[] vectorStrings = som.getStringValue("vector").split(",");			
				double[] vectorDoubles = new double[vectorStrings.length]; 

				for(int i = 0; i< vectorStrings.length; i++){
					vectorDoubles[i] = Double.parseDouble(vectorStrings[i]);
				}			

				// create data item (cluster)
				SOMDataItem dataItem = new SOMDataItem(color, values.toArray(new String[values.size()]), vectorDoubles);

				int x = som.getIntValue("x");
				int y = som.getIntValue("y");

				// array indexes start from zero
				dataset.addValue(x-1, y-1, dataItem);
			}
		}

		JFreeChart chart = BioChartFactory.createSOMChart(
				data.getName() + " SOM Chart",     // chart title
				dataset,                  // data
				true,                     // tooltips?
				false                     // URLs?
		);

		ChartPanel chartPanel = makePanel(chart);
		chartPanel.addChartMouseListener((SOMPlot)chart.getPlot());
		chartPanel.setFocusable(true);
		return chartPanel;
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return isTabular(bean) && bean.hasTypeTag(MicroarrayModule.TypeTags.SOM_CLUSTERED_EXPRESSION_VALUES);
	}
}
