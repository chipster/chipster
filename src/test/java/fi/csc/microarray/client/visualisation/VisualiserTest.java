package fi.csc.microarray.client.visualisation;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.WindowConstants;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.TestConstants;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
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
		this.manager = new DataManager();
		Session.getSession().setDataManager(this.manager);
		//testExpressionProfile fails when there is no application, but otherwise everything fails because of the missing configuration
		//Session.getSession().setClientApplication(new ClientContextUtil.SkeletonApplication());
		new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule").plugAll(manager, null);
	}

	public static void main(String[] args) throws IOException, Exception {
		new VisualiserTest().testHC();
	}
	
	@Test
	public void testSom() throws Exception {
		DataBean dataset = manager.createDataBean("SOM", getInputStream(TestConstants.SOM_CLUSTERED_RESOURCE));

		JComponent component = doVisualisation(MicroarrayModule.VisualisationMethods.SOM, dataset);
		makeFrame(component);
	}
	
	@Test
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
			DataBean tree = manager.createDataBean("HC clusters", getInputStream(trees[i]));
			DataBean heatmap = manager.createDataBean("heatmap", getInputStream(heatmaps[i]));
			tree.addLink(Link.DERIVATION, heatmap);

			JComponent component = MicroarrayModule.VisualisationMethods.HIERARCHICAL.getHeadlessVisualiser().getVisualisation(tree);
			makeFrame(component);
		}
	}
	
	@Test
	public void testVisualisations() throws Exception {

		String[] resources = new String[] {TestConstants.CDNA_RESOURCE, TestConstants.RESULSET_RESOURCE, TestConstants.AFFY_RESOURCE};
		for (String resource : resources) {
			try {
				DataBean bean = manager.createDataBean(resource, getInputStream(resource));
				
				doVisualisation(MicroarrayModule.VisualisationMethods.ARRAY_LAYOUT, bean);
				doVisualisation(BasicModule.VisualisationMethods.SPREADSHEET, bean);
				
			} catch (Exception e) {
				System.err.println("exception when processing " + resource + ": " + e.getMessage());
				throw e;
			}
		}
	}
	
	/**
	 * @throws Exception
	 */
	@Test
	public void testScatterplot() throws Exception {
		for (String resource : new String[] {TestConstants.FOUR_CHIPS_RESOURCE, TestConstants.SCATTER_HARDCASE1, TestConstants.SCATTER_HARDCASE2}) {
			try {
				DataBean dataBean = manager.createDataBean(resource, getInputStream(resource));
				
				Visualisation vis = MicroarrayModule.VisualisationMethods.SCATTERPLOT.getHeadlessVisualiser();
				
				vis.initialise(new InternalVisualisationFrame(FrameType.MAIN));
		
				vis.getVariablesFor(dataBean);
				vis.getParameterPanel();				
				JComponent component = vis.getVisualisation(dataBean);
				makeFrame(component);
				
			} catch (Exception e) {
				System.err.println("exception when processing " + resource + ": " + e.getMessage());
				throw e;
			}
		}
	}

	@Test
	public void testHistogram() throws Exception {
		
		DataBean dataBean = manager.createDataBean("Hist. data", getInputStream(TestConstants.FOUR_CHIPS_RESOURCE));		
		Variable[] variables = MicroarrayModule.VisualisationMethods.HISTOGRAM.getHeadlessVisualiser().getVariablesFor(dataBean);
		Assert.assertEquals(4, variables.length);
		Visualisation visualiser = MicroarrayModule.VisualisationMethods.HISTOGRAM.getHeadlessVisualiser();
		visualiser.getParameterPanel();		
		JComponent visualisation = visualiser.getVisualisation(dataBean);
		makeFrame(visualisation);
	}
	
	private FileInputStream getInputStream(String path) throws FileNotFoundException {
		//this.getClass().getResourceAsStream(path);
		File file = new File(path);
		return new FileInputStream(file);
	}

	@Test
	public void testExpressionProfile() throws Exception {
		
		DataBean dataBean = manager.createDataBean("Profiledata", getInputStream(TestConstants.FOUR_CHIPS_RESOURCE));
		JComponent visualisation = MicroarrayModule.VisualisationMethods.EXPRESSION_PROFILE.getHeadlessVisualiser().getVisualisation(dataBean);
		makeFrame(visualisation);
	}

	@Test
	public void testClusteredProfiles() throws Exception {
		//testfile missing
//		DataBean dataBean = manager.createDataBean("Profiledata", getInputStream(TestConstants.CLUSTERED_PROFILES_RESOURCE));
//		JComponent visualisation = MicroarrayModule.VisualisationMethods.CLUSTERED_PROFILES.getHeadlessVisualiser().getVisualisation(dataBean);
//		makeFrame(visualisation);
	}

	@Test
	public void testApplicabilityChecks() throws MicroarrayException, IOException {
		DataBean affyMicroarray = manager.createDataBean("", getInputStream(TestConstants.AFFY_RESOURCE));
		
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
