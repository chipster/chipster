package fi.csc.microarray.client.visualisation.methods;

import java.awt.Color;
import java.awt.geom.Rectangle2D;
import java.util.Collections;
import java.util.LinkedList;

import javax.swing.JComponent;

import org.apache.log4j.Logger;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.labels.StandardCategoryToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.DefaultCategoryItemRenderer;
import org.jfree.data.category.DefaultCategoryDataset;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.client.visualisation.TableAnnotationProvider;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel.SelectionChangeListener;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;

public class ExpressionProfile extends Visualisation implements SelectionChangeListener {
	
	private static final Logger logger = Logger.getLogger(ExpressionProfile.class);

	public ExpressionProfile(VisualisationFrame frame) {
		super(frame);
	}

	// color scale red->yellow->blue
	private static final float START_COLOR_R = 0.0f;
	private static final float START_COLOR_G = 0.0f;
	private static final float START_COLOR_B = 1.0f;
	private static final float MIDDLE_COLOR_R = 1.0f;
	private static final float MIDDLE_COLOR_G = 1.0f;
	private static final float MIDDLE_COLOR_B = 0.0f;
	private static final float END_COLOR_R = 1.0f;
	private static final float END_COLOR_G = 0.0f;
	private static final float END_COLOR_B = 0.0f;

	private static Color getColor(float position) {
		boolean beforeMiddle = position < 0.5f;
		float gradientPosition = beforeMiddle ? position * 2.0f : (position-0.5f) * 2.0f;
		float invGradientPosition = 1.0f - gradientPosition;

		float r, g, b;
		if (beforeMiddle) {
			r = START_COLOR_R * invGradientPosition + MIDDLE_COLOR_R * gradientPosition;
			g = START_COLOR_G * invGradientPosition + MIDDLE_COLOR_G * gradientPosition;
			b = START_COLOR_B * invGradientPosition + MIDDLE_COLOR_B * gradientPosition;
		} else {
			r = MIDDLE_COLOR_R * invGradientPosition + END_COLOR_R * gradientPosition;
			g = MIDDLE_COLOR_G * invGradientPosition + END_COLOR_G * gradientPosition;
			b = MIDDLE_COLOR_B * invGradientPosition + END_COLOR_B * gradientPosition;
		}

		return new Color(r, g, b);			
	}		
	
	public static class ProfileRow implements Comparable<ProfileRow> {
		int series;
		Float value;
		Color color;				

		public int compareTo(ProfileRow row) {
			return value.compareTo(row.value);
		}
	}
	
	/**
	 * The user can change the names of the columns in the phenodata and they aren't necessary 
	 * individual as they might show for example group information. To keep the values of separate chips separated
	 * the equals comparison has to be done according internal column name.
	 * 
	 * @author klemela
	 */
	public static class IndividualizedColumn implements Comparable<IndividualizedColumn> {
		private String viewName;
		private String internalName;
		public IndividualizedColumn(String internalName, String viewName){
			this.internalName = internalName;
			this.viewName = viewName;
		}
		
		public int compareTo(IndividualizedColumn o) {			
			return internalName.compareTo(o.internalName);
		}				
		
		@Override
		public String toString(){			
			return viewName;
		}
		
		@Override
		public boolean equals(Object o){
			if (o instanceof IndividualizedColumn) {
				IndividualizedColumn ic = (IndividualizedColumn) o;
				
				return internalName.equals(ic.internalName);
			} else {
				return false;
			}
		}
		
		@Override
		public int hashCode(){
			return internalName.hashCode();
		}
	}
	
	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		
		TableAnnotationProvider annotationProvider = new TableAnnotationProvider(data);

		// fetch data
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		
		Table samples = data.queryFeatures("/column/*").asTable();
		
		// read through data
		LinkedList<ProfileRow> rows = new LinkedList<ProfileRow>();
		int rowNumber = 0;
		while (samples.nextRow()) {
			boolean firstSample = true;
			for (String sample : samples.getColumnNames()) {
				// collect all chip columns
				if (sample.startsWith("chip.")) {
					
					// order by first chip
					if (firstSample) {
						ProfileRow row = new ProfileRow();
						row.value = samples.getFloatValue(sample);
						row.series = rowNumber;
						firstSample = false;
						rows.add(row);
					}

					// insert into dataset 
					String sampleName = data.queryFeatures("/phenodata/linked/describe/" + sample.substring("chip.".length())).asString();
					String rowName = annotationProvider.getAnnotatedRowname(samples.getStringValue(" "));
					IndividualizedColumn column = new IndividualizedColumn(sample, sampleName);
					dataset.addValue((double)samples.getFloatValue(sample), rowName, column);
				}
			}
			rowNumber++;
		}

		JFreeChart chart = createProfileChart(dataset, rows, data.getName());
		//return makeSelectablePanel(chart, this);
		return makePanel(chart);
	}

	public static JFreeChart createProfileChart(DefaultCategoryDataset dataset, LinkedList<ProfileRow> rows, String name) {
		// create renderer
        DefaultCategoryItemRenderer renderer = new DefaultCategoryItemRenderer();        
        renderer.setToolTipGenerator(new StandardCategoryToolTipGenerator());
        renderer.setShapesVisible(false);
        
        // generate colors
        Collections.sort(rows);
		float position = 0.0f;
		float step = 1.0f / ((float)rows.size());
		for (ProfileRow row : rows) {
			row.color = getColor(position);
			renderer.setSeriesPaint(row.series, getColor(position));
			position += step;
		}

		// draw plot
        CategoryAxis categoryAxis = new CategoryAxis("sample");
        categoryAxis.setCategoryLabelPositions(CategoryLabelPositions.UP_45);
        categoryAxis.setUpperMargin(0.0);
        categoryAxis.setLowerMargin(0.0);
        ValueAxis valueAxis = new NumberAxis("expression");
        CategoryPlot plot = new CategoryPlot(dataset, categoryAxis, valueAxis, renderer);
        plot.setOrientation(PlotOrientation.VERTICAL);
        
        JFreeChart chart = new JFreeChart("Expression profile for " + name, 
        		JFreeChart.DEFAULT_TITLE_FONT, plot, false);
        
		return chart;
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		boolean isTabular = VisualisationMethod.SPREADSHEET.getHeadlessVisualiser().canVisualise(bean);
		if (isTabular) {
			Table chips = bean.queryFeatures("/column/chip.*").asTable();
			return chips != null && chips.getColumnNames().length > 1;
		}
		return false;
	}

	public void selectionChanged(Rectangle2D.Double selection) {
		if(selection == null){
			logger.debug("Empty selection");
		} else if (selection.getWidth() == 0 && selection.getHeight() == 0){
			logger.debug("Single click at (" + selection.getX() + ", " + selection.getY());
		} else {
			logger.debug("Area selection: " + 
					selection.getMinX() + " < X < " + selection.getMaxX() + " and " + 
					selection.getMinY() + " < Y < " + selection.getMaxY());
		}
	}	
}
