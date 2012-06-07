package fi.csc.microarray.client.visualisation;

import java.io.IOException;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.TestConstants;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class VisualiserTest {

	private DataManager manager;
	
	public VisualiserTest() throws Exception {
		DirectoryLayout.uninitialise();
		DirectoryLayout.initialiseSimpleLayout().getConfiguration();			
		this.manager = new DataManager();
		new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule").plugAll(manager, null);
	}

	public static void main(String[] args) throws IOException, Exception {
		new VisualiserTest().testHC();
	}
	
	@Test(groups = {"smoke"} )
	public void testSom() throws Exception {
		DataBean dataset = manager.createDataBean("SOM", this.getClass().getResourceAsStream(TestConstants.SOM_CLUSTERED_RESOURCE));

		JComponent component = doVisualisation(MicroarrayModule.VisualisationMethods.SOM, dataset);
		makeFrame(component);
	}
	
	@Test(groups = {"smoke"} )
	public void testHC() throws Exception {
		String[] trees = new String[] {
				TestConstants.HIERARCHICAL_CLUSTERED_RESOURCE,
				TestConstants.HIERARCHICAL_CLUSTERED_ILLUMINA_RESOURCE,
				TestConstants.HIERARCHICAL_CLUSTERED_AGILENT_RESOURCE
		};
		String[] heatmaps = new String[] {
				TestConstants.HIERARCHICAL_CLUSTERED_HEATMAP_RESOURCE,
				TestConstants.HIERARCHICAL_CLUSTERED_ILLUMINA_HEATMAP_RESOURCE,
				TestConstants.HIERARCHICAL_CLUSTERED_AGILENT_HEATMAP_RESOURCE
		};
		
		for (int i = 0; i < trees.length; i++) {
			DataBean tree = manager.createDataBean("HC clusters", this.getClass().getResourceAsStream(trees[i]));
			DataBean heatmap = manager.createDataBean("heatmap", this.getClass().getResourceAsStream(heatmaps[i]));
			tree.addLink(Link.DERIVATION, heatmap);

			JComponent component = MicroarrayModule.VisualisationMethods.HIERARCHICAL.getHeadlessVisualiser().getVisualisation(tree);
			makeFrame(component);
		}
	}
	
	@Test(groups = {"smoke"} )
	public void testVisualisations() throws Exception {

		String[] resources = new String[] {TestConstants.CDNA_RESOURCE, TestConstants.RESULSET_RESOURCE, TestConstants.AFFY_RESOURCE};
		for (String resource : resources) {
			try {
				DataBean bean = manager.createDataBean(resource, this.getClass().getResourceAsStream(resource));
				
				doVisualisation(MicroarrayModule.VisualisationMethods.ARRAY_LAYOUT, bean);
				doVisualisation(BasicModule.VisualisationMethods.SPREADSHEET, bean);
				
			} catch (Exception e) {
				System.err.println("exception when processing " + resource + ": " + e.getMessage());
				throw e;
			}
		}
	}
	
	@Test(groups = {"smoke"} )
	public void testScatterplot() throws Exception {
		for (String resource : new String[] {TestConstants.FOUR_CHIPS_RESOURCE, TestConstants.SCATTER_HARDCASE1, TestConstants.SCATTER_HARDCASE2}) {
			try {
				DataBean dataBean = manager.createDataBean(resource, this.getClass().getResourceAsStream(resource));
				
				MicroarrayModule.VisualisationMethods.SCATTERPLOT.getHeadlessVisualiser().getVariablesFor(dataBean);
				MicroarrayModule.VisualisationMethods.SCATTERPLOT.getHeadlessVisualiser().getParameterPanel();
				JComponent visualisation = MicroarrayModule.VisualisationMethods.SCATTERPLOT.getHeadlessVisualiser().getVisualisation(dataBean);
				makeFrame(visualisation);
				
			} catch (Exception e) {
				System.err.println("exception when processing " + resource + ": " + e.getMessage());
				throw e;
			}
		}
	}

	@Test(groups = {"smoke"} )
	public void testHistogram() throws Exception {
		
		DataBean dataBean = manager.createDataBean("Hist. data", this.getClass().getResourceAsStream(TestConstants.FOUR_CHIPS_RESOURCE));		
		Variable[] variables = MicroarrayModule.VisualisationMethods.HISTOGRAM.getHeadlessVisualiser().getVariablesFor(dataBean);
		Assert.assertEquals(variables.length, 4);
		Visualisation visualiser = MicroarrayModule.VisualisationMethods.HISTOGRAM.getHeadlessVisualiser();
		visualiser.getParameterPanel();		
		JComponent visualisation = visualiser.getVisualisation(dataBean);
		makeFrame(visualisation);
	}
	
	@Test(groups = {"smoke"} )
	public void testExpressionProfile() throws Exception {
		
		DataBean dataBean = manager.createDataBean("Profiledata", this.getClass().getResourceAsStream(TestConstants.FOUR_CHIPS_RESOURCE));
		JComponent visualisation = MicroarrayModule.VisualisationMethods.EXPRESSION_PROFILE.getHeadlessVisualiser().getVisualisation(dataBean);
		makeFrame(visualisation);
	}

	@Test(groups = {"smoke"} )
	public void testClusteredProfiles() throws Exception {
		
		DataBean dataBean = manager.createDataBean("Profiledata", this.getClass().getResourceAsStream(TestConstants.CLUSTERED_PROFILES_RESOURCE));
		JComponent visualisation = MicroarrayModule.VisualisationMethods.CLUSTERED_PROFILES.getHeadlessVisualiser().getVisualisation(dataBean);
		makeFrame(visualisation);
	}

	@Test(groups = {"smoke"} )
	public void testApplicabilityChecks() throws MicroarrayException, IOException {
		DataBean affyMicroarray = manager.createDataBean("", this.getClass().getResourceAsStream(TestConstants.AFFY_RESOURCE));
		
		for (VisualisationMethod method : Session.getSession().getVisualisations().getVisualisationMethods()) {
			if (method.isApplicableTo(affyMicroarray)) {
				if (method == MicroarrayModule.VisualisationMethods.SOM || method == MicroarrayModule.VisualisationMethods.HIERARCHICAL
						|| method == MicroarrayModule.VisualisationMethods.EXPRESSION_PROFILE || method == BasicModule.VisualisationMethods.SHOW_IMAGE
						|| method == BasicModule.VisualisationMethods.WEBVIEW || method == BasicModule.VisualisationMethods.VIEW_TEXT) {
					
					Assert.fail("method " + method.getName() + " should not be applicable to " + affyMicroarray.getName());
				}
			}
		}
	}
	
	private JComponent doVisualisation(VisualisationMethod method, DataBean bean) throws Exception {
		method.getHeadlessVisualiser().getParameterPanel();
		return method.getHeadlessVisualiser().getVisualisation(bean);
	}


	private void makeFrame(JComponent visualisation) {
		JFrame frame = new JFrame();
		frame.add(visualisation);
		frame.pack();
		frame.setVisible(true);
		frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE);
	}
}
