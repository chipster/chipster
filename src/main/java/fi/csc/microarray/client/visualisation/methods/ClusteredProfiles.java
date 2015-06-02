package fi.csc.microarray.client.visualisation.methods;

import java.awt.GridLayout;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JPanel;

import org.jfree.chart.JFreeChart;
import org.jfree.chart.MultiChartPanel;
import org.jfree.data.category.DefaultCategoryDataset;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.TableAnnotationProvider;
import fi.csc.microarray.client.visualisation.Visualisation;
import fi.csc.microarray.client.visualisation.VisualisationFrame;
import fi.csc.microarray.client.visualisation.methods.ExpressionProfile.ProfileRow;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class ClusteredProfiles extends Visualisation {

	public void initialise(VisualisationFrame frame) throws Exception {
		super.initialise(frame);
	}

	@Override
	public JComponent getVisualisation(DataBean data) throws Exception {
		
		// create a local copy of the phenodata to speed up initialization over
		// slower connections (takes several minutes otherwise)
		DataBean metadata = LinkUtils.retrieveInherited(data, Link.ANNOTATION);
		Session.getSession().getDataManager().getLocalFile(metadata);
		
		// count clusters
		int clusterCount = 0;
		Iterable<String> clusters = data.queryFeatures("/column/cluster").asStrings();
		for (String cluster : clusters) {
			int num = Integer.parseInt(cluster);
			if (num > clusterCount) {
				clusterCount = num;
			}
		}
		
		// initialise cluster collector
	    ArrayList<JFreeChart> charts = new ArrayList<JFreeChart>(clusterCount);

		// fetch datas
		LinkedList<DefaultCategoryDataset> datasets = new LinkedList<DefaultCategoryDataset>();
		LinkedList<LinkedList<ProfileRow>> rows = new LinkedList<LinkedList<ProfileRow>>();
		for (int i = 0; i < clusterCount; i++) {
			datasets.add(new DefaultCategoryDataset());
			rows.add(new LinkedList<ProfileRow>());
		}
		
		
		TableAnnotationProvider annotationProvider = new TableAnnotationProvider(data);
		
		try (Table samples = data.queryFeatures("/column/*").asTable()) {		
			int[] rowNumbers = new int[clusterCount];
			while (samples.nextRow()) {
				int cluster = samples.getIntValue("cluster");
				boolean firstSample = true;
				for (String sample : samples.getColumnNames()) {
					if (sample.startsWith("chip.")) {

						// order by first chip
						if (firstSample) {
							ProfileRow row = new ProfileRow();
							row.value = samples.getFloatValue(sample);
							row.series = rowNumbers[cluster-1];
							firstSample = false;
							rows.get(cluster-1).add(row);
						}
						String sampleName = data.queryFeatures("/phenodata/linked/describe/" + sample.substring("chip.".length())).asString();
						datasets.get(cluster-1).addValue((double)samples.getFloatValue(sample), samples.getStringValue(" "), sampleName);
						String rowName = annotationProvider.getAnnotatedRowname(samples.getStringValue(" "));
						datasets.get(cluster-1).addValue((double)samples.getFloatValue(sample), rowName, sampleName);
					}
				}
				rowNumbers[cluster-1]++;
			}
		}
		
		// draw plots
		for (int i = 0; i < clusterCount; i++) {
			//TODO createProfileChart should allow static use
			// FIXME it really should
			ExpressionProfile ep = new ExpressionProfile();
			ep.initialise(null);
			charts.add(ep.createProfileChart(datasets.get(i), rows.get(i), "Cluster " + (i+1)));
		}
		
		return this.makePanel(charts);
	}
	
	protected JComponent makePanel(List<JFreeChart> charts) {
		
		// calculate coordinates
		int slots = (int)Math.ceil(Math.sqrt(charts.size()));
		int height = slots;
		int width = slots*(slots-1) >= charts.size() ? slots-1 : slots; // shave off from width if possible

		// place panels
		JPanel panel = new JPanel(new GridLayout(width, height));
        LinkedList<MultiChartPanel> panels = new LinkedList<MultiChartPanel>();
		for (JFreeChart chart: charts) {
			MultiChartPanel multiPanel = new MultiChartPanel(chart, panels);
			panels.add(multiPanel);
			panel.add(multiPanel);
		}
		//panel.setPreferredSize(size);
		
		return panel;
	}

	@Override
	public boolean canVisualise(DataBean bean) throws MicroarrayException {
		return isTabular(bean) && bean.hasTypeTag(MicroarrayModule.TypeTags.NORMALISED_EXPRESSION_VALUES) && bean.hasTypeTag(MicroarrayModule.TypeTags.CLUSTERED_EXPRESSION_VALUES);
	}
	
}
