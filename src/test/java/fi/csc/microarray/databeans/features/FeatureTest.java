package fi.csc.microarray.databeans.features;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.TestConstants;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.fs.FSDataManager;
import fi.csc.microarray.module.DefaultModules;

public class FeatureTest {

	private DataManager manager;

	public FeatureTest() throws IOException, IllegalConfigurationException {
		DirectoryLayout.initialiseUnitTestLayout();
		this.manager = new FSDataManager();
		DefaultModules.getDefaultModules().plugFeatures(manager);
	}

	@Test(groups = {"unit"} )
	public void testEmbeddedBinary() throws IOException, MicroarrayException {
		DataBean affyMicroarray = manager.createDataBean("affy.cel", new FileInputStream(TestConstants.AFFY_RESOURCE));
		DataBean binAffyMicroarray = manager.createDataBean("bin_affy.cel", new FileInputStream(TestConstants.BIN_AFFY_RESOURCE));
		
		Assert.assertFalse(affyMicroarray.queryFeatures("/embedded-binary-content/").exists());
		Assert.assertTrue(binAffyMicroarray.queryFeatures("/embedded-binary-content/").exists());
	}

	@Test(groups = {"unit"} )
	public void testPhenodataFeatures() throws IOException, MicroarrayException {
		DataBean data = manager.createDataBean("filtered.tsv", new FileInputStream(TestConstants.FOUR_CHIPS_RESOURCE));
		DataBean phenoData= manager.createDataBean("phenodata.tsv", new FileInputStream(TestConstants.FOUR_CHIPS_PHENODATA_RESOURCE));
		phenoData.addLink(Link.ANNOTATION, data);
		Assert.assertTrue(phenoData.queryFeatures("/phenodata/").exists());
		Assert.assertTrue(phenoData.queryFeatures("/phenodata/is_complete").exists());
		Assert.assertEquals(phenoData.queryFeatures("/phenodata/describe/microarray1.cel").asString(), "GSM11814.cel");
		Assert.assertEquals(phenoData.queryFeatures("/phenodata/describe/not_in_phenodata").asString(), "not_in_phenodata");
		
		Assert.assertTrue(data.queryFeatures("/phenodata/linked/").exists());
		Assert.assertTrue(data.queryFeatures("/phenodata/linked/is_complete").exists());
	}

	@Test(groups = {"unit"} )
	public void testModifiers() throws IOException, MicroarrayException {
		DataBean affyMicroarray = manager.createDataBean("affy.cel", new FileInputStream(TestConstants.AFFY_RESOURCE));
		QueryResult feature = affyMicroarray.queryFeatures("log(/normalised-expression)");
		Assert.assertTrue(feature.exists());

		// check that can be iterated twice
		float last = 1;
		for (float f : feature.asFloats()) {
			last = f;
		}
		for (float f : feature.asFloats()) {
			last = f;
		}
		Assert.assertEquals(last, 7.426265f);
		
		QueryResult doubleFeature = affyMicroarray.queryFeatures("log(log(/normalised-expression))");
		Assert.assertTrue(doubleFeature.exists());

	}
	
	public static void main(String[] args) throws IOException, MicroarrayException, IllegalConfigurationException {
		new FeatureTest().testRowCount();
	}

	@Test(groups = {"unit"} )
	public void testRowCount() throws MicroarrayException, FileNotFoundException {
		DataBean affyMicroarray = manager.createDataBean("affy.cel", new FileInputStream(TestConstants.AFFY_RESOURCE));
		Assert.assertEquals(affyMicroarray.queryFeatures("/rowcount/max/10").asFloat(), 10f);
		Assert.assertEquals(affyMicroarray.queryFeatures("/rowcount/max/1000000").asFloat(), 15876f);
	}
	
	@Test(groups = {"unit"} )
	public void testTableColumnIterable() throws MicroarrayException, FileNotFoundException {
		DataBean affyMicroarray = manager.createDataBean("affy.cel", new FileInputStream(TestConstants.AFFY_RESOURCE));
		QueryResult mean = affyMicroarray.queryFeatures("/column/MEAN");
		Iterable[] iterables = new Iterable[] {
				mean.asFloats(),
				mean.asStrings()
		};
		for (Iterable iterable : iterables) {
			Iterator iterator = iterable.iterator();
			
			// can we call next() directly
			iterator.next(); 
			
			// can we iterate over rest of it
			while (iterator.hasNext()) {
				iterator.next();
			}
		}
	}
	
	@Test(groups = {"unit"} )
	public void testFeatures() throws MicroarrayException, IOException {
		DataBean affyMicroarray = manager.createDataBean("affy.cel", new FileInputStream(TestConstants.AFFY_RESOURCE));
		for (String feature : new String [] {"/normalised-expression", "/column/MEAN"}) {			
			Assert.assertNotNull(affyMicroarray.queryFeatures(feature).asFloats(),"error in " + feature);
			for (float f : affyMicroarray.queryFeatures(feature).asFloats()) {
				Assert.assertTrue(f > 0.0, "illegal value: " + f + " in " + feature);
			}
		}
		Assert.assertTrue(affyMicroarray.queryFeatures("/column/MEAN").asFloats().iterator().next() == 190f);
		float last = 0;
		for (float f : affyMicroarray.queryFeatures("/column/MEAN").asFloats()) {
			last = f;
		}
		Assert.assertTrue(last == 172f);
		Assert.assertTrue(VisualisationMethod.ARRAY_LAYOUT.isApplicableTo(affyMicroarray));

		// SOM
		DataBean somData = manager.createDataBean("som.tsv", new FileInputStream(TestConstants.SOM_CLUSTERED_RESOURCE));

		Assert.assertTrue(VisualisationMethod.SOM.isApplicableTo(somData));
		Table som = somData.queryFeatures("/clusters/som").asTable();
		Assert.assertNotNull(som);
		Assert.assertEquals(som.getColumnCount(), 5); 
		
		// hierarchical clustering
		DataBean hcTree = manager.createDataBean("hs.tre", new FileInputStream(TestConstants.HIERARCHICAL_CLUSTERED_RESOURCE));
		
		DataBean hcHeatmap = manager.createDataBean("hc.tsv", new FileInputStream(TestConstants.HIERARCHICAL_CLUSTERED_HEATMAP_RESOURCE));
		

		hcTree.addLink(Link.DERIVATION, hcHeatmap);

		Assert.assertTrue(VisualisationMethod.HIERARCHICAL.isApplicableTo(hcTree));
		Table heatmap = hcTree.queryFeatures("/clusters/hierarchical/heatmap").asTable();
		Assert.assertNotNull(heatmap);
		heatmap.nextRow();
		Assert.assertNotNull(heatmap.getValue(" "));
		
		String tree = hcTree.queryFeatures("/clusters/hierarchical/tree").asStrings().iterator().next();
		Assert.assertNotNull(tree);
	}
}
