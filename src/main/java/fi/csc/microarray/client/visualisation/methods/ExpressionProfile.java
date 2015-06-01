package fi.csc.microarray.client.visualisation.methods;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.geom.Line2D;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Set;

import javax.swing.JComponent;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.CategoryAxis;
import org.jfree.chart.axis.CategoryLabelPositions;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.labels.StandardCategoryToolTipGenerator;
import org.jfree.chart.plot.CategoryPlot;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.renderer.category.CategoryItemRenderer;
import org.jfree.chart.renderer.category.DefaultCategoryItemRenderer;
import org.jfree.data.category.CategoryDataset;
import org.jfree.data.category.DefaultCategoryDataset;

import fi.csc.microarray.client.selection.SelectionEvent;
import fi.csc.microarray.client.selection.IntegratedSelectionManager;
import fi.csc.microarray.client.visualisation.SelectionList;
import fi.csc.microarray.client.visualisation.TableAnnotationProvider;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.SelectableChartPanel.SelectionChangeListener;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;

public class ExpressionProfile extends ChipVisualisation 
implements PropertyChangeListener, SelectionChangeListener {
	
	//private static final Logger logger = Logger.getLogger(ExpressionProfile.class);
	
	private LinkedList<ProfileRow> rows;
	private SelectableChartPanel selectableChartPanel;
	private  JPanel paramPanel;
	private SelectionList list;

	private CategoryPlot plot;
	
	private DataBean data;
	
	//selection indexes in order of the original data
	private Set<Integer> selectedIndexes = new HashSet<Integer>();

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}
	
	@Override
	public JPanel getParameterPanel() {
		if (paramPanel == null) {
			paramPanel = new JPanel();
			paramPanel.setPreferredSize(Visualisation.PARAMETER_SIZE);
			paramPanel.setLayout(new BorderLayout());

			list = new SelectionList();

			JTabbedPane tabPane = new JTabbedPane();
			tabPane.addTab("Selected", list);

			paramPanel.add(tabPane, BorderLayout.CENTER);
		}
		return paramPanel;
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

	/**
	 * Calculates Color that corresponds to position in the gradient.
	 */
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

		return new Color(limit(r), limit(g), limit(b));			
	}		

	/**
	 * Safe guards float values between 0-1.
	 */
	private static float limit(float f) {
		return Math.max(0.0f,  Math.min(1.0f, f));
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
		
		this.data = data;
		
		application.addClientEventListener(this);
		
		JFreeChart chart = createProfileChart(createDataset(), rows,  data.getName());
		
		updateSelectionsFromApplication(false);
		
		selectableChartPanel = new SelectableChartPanel(chart, this); 
		return selectableChartPanel;
	}
	
	private CategoryDataset createDataset() throws MicroarrayException{
		TableAnnotationProvider annotationProvider = new TableAnnotationProvider(data);

		
		// make a lookup table for sample names
		HashMap<String, String> sampleNameLookup = new HashMap<String, String>();
		try (Table chips = data.queryFeatures("/column/chip.*").asTable()) {
			for (String chip : chips.getColumnNames()) {
				String sampleName = data.queryFeatures("/phenodata/linked/describe/" + chip.substring("chip.".length())).asString();
				sampleNameLookup.put(chip, sampleName);
			}
		}
		
		// fetch data
		DefaultCategoryDataset dataset = new DefaultCategoryDataset();
		try (Table samples = data.queryFeatures("/column/*").asTable()) {

			// read through data
			rows = new LinkedList<ProfileRow>();
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
						String sampleName = sampleNameLookup.get(sample);
						String rowName = annotationProvider.getAnnotatedRowname(samples.getStringValue(" "));
						IndividualizedColumn column = new IndividualizedColumn(sample, sampleName);
						dataset.addValue((double)samples.getFloatValue(sample), rowName, column);
					}
				}
				rowNumber++;
			}
		}

		return dataset;
	}
	
	private CategoryItemRenderer createRenderer(List<ProfileRow> rows){
		// create renderer
        DefaultCategoryItemRenderer renderer = new DefaultCategoryItemRenderer();        
        renderer.setToolTipGenerator(new StandardCategoryToolTipGenerator());
        renderer.setShapesVisible(false);
        
        // generate colors
        Collections.sort(rows);
		float position = 0.0f;
		float step = 1.0f / ((float)rows.size());
			
		for (ProfileRow row : rows) {
			if(selectedIndexes.contains(row.series)){
				renderer.setSeriesPaint(row.series, Color.black);
			} else {
				row.color = getColor(position);
				renderer.setSeriesPaint(row.series, getColor(position));
			}
			position += step;
		}
		
		//List isn't initialised, if visualising clustered profiles
		if(list != null){
			list.setSelectedRows(selectedIndexes, this, false, data);
		}
		
		return renderer;
	}

	public JFreeChart createProfileChart(CategoryDataset categoryDataset, List<ProfileRow> rows, String name) throws MicroarrayException {
						
		// draw plot
        CategoryAxis categoryAxis = new CategoryAxis("sample");
        categoryAxis.setCategoryLabelPositions(CategoryLabelPositions.UP_45);
        categoryAxis.setUpperMargin(0.0);
        categoryAxis.setLowerMargin(0.0);
        ValueAxis valueAxis = new NumberAxis("expression");
        plot = new CategoryPlot(categoryDataset, categoryAxis, valueAxis, createRenderer(rows));
        plot.setOrientation(PlotOrientation.VERTICAL);
        
        JFreeChart chart = new JFreeChart("Expression profile for " + name, 
        		JFreeChart.DEFAULT_TITLE_FONT, plot, false);
        
		return chart;
	}

	public void selectionChanged(Rectangle2D.Double selection) {
		
		if(selection == null){
			selectedIndexes.clear();
		} else {
						
			//Goes through all the lines in profile, some optimisation can be done if 
			//this is too slow
			
			Set<Integer>  newSelection = new HashSet<Integer>();

			CategoryDataset dataset = plot.getDataset();
			List colKeys = dataset.getColumnKeys();
			List rowKeys = dataset.getRowKeys();

			for(int x = 0; x < dataset.getColumnCount() - 1; x++){
				for (int y = 0; y < dataset.getRowCount(); y++){

					Point2D start = new Point2D.Double(x, 
							(Double) dataset.getValue((Comparable)rowKeys.get(y), 
									(Comparable) colKeys.get(x)));
					
					Point2D end = new Point2D.Double(x + 1, 
							(Double) dataset.getValue((Comparable)rowKeys.get(y), 
									(Comparable) colKeys.get(x + 1)));
					
					Line2D.Double line = new Line2D.Double(start, end);

					if(selection.intersectsLine(line)){
						newSelection.add(y);
					}
				}								
			}
			
			//New selections can't be put directly to selectedIndexes, because every other occurrence
			//of line inside selection rectangle would undo selection
			
			for(Integer row: newSelection){
				if(selectedIndexes.contains(row)){
					selectedIndexes.remove(row);
				} else {
					selectedIndexes.add(row);
				}
			}
		}
		
		this.list.setSelectedRows(selectedIndexes, this, true, data);
		
		plot.setRenderer(createRenderer(rows));
	}

	public void propertyChange(PropertyChangeEvent evt) {
		if (evt instanceof SelectionEvent && evt.getSource() != this && ((SelectionEvent) evt).getData() == data) {

			updateSelectionsFromApplication(false);
		}
	}	
	
	protected void updateSelectionsFromApplication(boolean dispatchEvent) {
		IntegratedSelectionManager manager = application.getSelectionManager().getSelectionManager(data);

		selectedIndexes.clear();
		for (int i : manager.getSelectionAsRows()){
			selectedIndexes.add(i);			
		}

		list.setSelectedRows(selectedIndexes, this, dispatchEvent, data);
		
		plot.setRenderer(createRenderer(rows));
	}
}
