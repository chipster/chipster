package fi.csc.microarray.client.visualisation.methods;

import javax.swing.JComponent;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.biobeans.BioBean;
import fi.csc.microarray.proto.MyPlot;
import fi.csc.microarray.util.FloatArrayList;

public class ArrayLayout extends Visualisation {

public ArrayLayout(VisualisationFrame frame) {
		super(frame);
	}

	//	@SuppressWarnings("all") // deprecated JFreeChart API in use...
	@SuppressWarnings("deprecation")
	@Override
	public JComponent getVisualisation(DataBean dataBean) throws Exception {
		
		// get intensities
		FloatArrayList transformedIntensities = new FloatArrayList();
		
		// use log to tune down outliers
		for (float intensity : dataBean.queryFeatures("log(/normalised-expression)").asFloats()) {
			
			if (intensity <= 0.0) {
				intensity = Float.MIN_VALUE;
			}
			transformedIntensities.add(intensity); 
		}
		
		// scale to 0...1
		float max = transformedIntensities.max();
		
		for (int c = 0; c < transformedIntensities.size(); c++) {
			transformedIntensities.set(c, transformedIntensities.get(c) / max);
		}
		
		// build heatmap
		double[] xs = new double[transformedIntensities.size()];
		double[] ys = new double[transformedIntensities.size()];
		double[] zs = new double[transformedIntensities.size()];
		int i = 0;
		int[] dimensions = inferDimensions(transformedIntensities);
		for (int x = 0; x < dimensions[0]; x++) {
			for (int y = 0; y < dimensions[1]; y++) {
				xs[i] = x;
				ys[i] = y;
				zs[i] = transformedIntensities.getFloat((dimensions[1]-y-1)*dimensions[0]+x); // invert Y for matching coordinate systems 
				i++;				
			}
		}

		org.jfree.data.contour.ContourDataset data = new org.jfree.data.contour.DefaultContourDataset("data",
				org.jfree.data.contour.DefaultContourDataset.formObjectArray(xs), 
				org.jfree.data.contour.DefaultContourDataset.formObjectArray(ys), 
				org.jfree.data.contour.DefaultContourDataset.formObjectArray(zs));
		org.jfree.chart.axis.ColorBar colors = new org.jfree.chart.axis.ColorBar("Expression");
		colors.setColorPalette(new org.jfree.chart.plot.GreyPalette());
		colors.autoAdjustRange();

		ValueAxis xAxis = new NumberAxis("x");
		ValueAxis yAxis = new NumberAxis("y");

		MyPlot plot = new MyPlot(data, xAxis, yAxis, colors);
		JFreeChart chart = new JFreeChart(plot);
		chart.setTitle(dataBean.getName());
		
		return makePanel(chart);
	}

	public static int[] inferDimensions(FloatArrayList intensities) throws MicroarrayException {
		// infer size
		int i = (int)Math.sqrt(intensities.size());
		// try different sizes - this should end at least when i = 1
		while (intensities.size() % i != 0) {
			i--;
		}

		return new int[] {i, intensities.size() / i};
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		boolean isTabular = VisualisationMethod.SPREADSHEET.isApplicableTo(bean);
		return isTabular && new BioBean(bean).getColorCount() == 1;
	}

}
