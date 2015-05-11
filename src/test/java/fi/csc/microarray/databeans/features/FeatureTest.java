package fi.csc.microarray.databeans.features;

import java.io.FileInputStream;
import java.io.IOException;
import java.util.Iterator;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.TestConstants;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.basic.BasicModule;
import fi.csc.microarray.module.chipster.MicroarrayModule;

public class FeatureTest {

	private DataManager manager;

	public FeatureTest() throws Exception {
		DirectoryLayout.uninitialise();
		DirectoryLayout.initialiseUnitTestLayout();
		this.manager = new DataManager();
		new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule").plugAll(manager, null);
	}

	@Test
	public void testPhenodataFeatures() throws IOException, MicroarrayException {
		DataBean data = manager.createDataBean("filtered.tsv", new FileInputStream(TestConstants.FOUR_CHIPS_RESOURCE));
		DataBean phenoData= manager.createDataBean("phenodata.tsv", new FileInputStream(TestConstants.FOUR_CHIPS_PHENODATA_RESOURCE));
		phenoData.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
		phenoData.addLink(Link.ANNOTATION, data);
		Assert.assertTrue(phenoData.queryFeatures("/phenodata/").exists());
		Assert.assertTrue(phenoData.queryFeatures("/phenodata/is_complete").exists());
		Assert.assertEquals(phenoData.queryFeatures("/phenodata/describe/microarray1.cel").asString(), "GSM11814.cel");
		Assert.assertEquals(phenoData.queryFeatures("/phenodata/describe/not_in_phenodata").asString(), "not_in_phenodata");
		
		Assert.assertTrue(data.queryFeatures("/phenodata/linked/").exists());
		Assert.assertTrue(data.queryFeatures("/phenodata/linked/is_complete").exists());
	}

	@Test
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
		Assert.assertEquals(last, 7.426265f, 0.01f);
		
		QueryResult doubleFeature = affyMicroarray.queryFeatures("log(log(/normalised-expression))");
		Assert.assertTrue(doubleFeature.exists());

	}
	
	@Test
	public void testTableColumnIterable() throws MicroarrayException, IOException {
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
	
	@Test
	public void testFeatures() throws MicroarrayException, IOException {
		DataBean affyMicroarray = manager.createDataBean("affy.cel", new FileInputStream(TestConstants.AFFY_RESOURCE));
		affyMicroarray.addTypeTag(MicroarrayModule.TypeTags.RAW_AFFYMETRIX_EXPRESSION_VALUES);
		affyMicroarray.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
		for (String feature : new String [] {"/normalised-expression", "/column/MEAN"}) {			
			Assert.assertNotNull("error in " + feature, affyMicroarray.queryFeatures(feature).asFloats());
			for (float f : affyMicroarray.queryFeatures(feature).asFloats()) {
				Assert.assertTrue("illegal value: " + f + " in " + feature, f > 0.0);
			}
		}
		Assert.assertTrue(affyMicroarray.queryFeatures("/column/MEAN").asFloats().iterator().next() == 190f);
		float last = 0;
		for (float f : affyMicroarray.queryFeatures("/column/MEAN").asFloats()) {
			last = f;
		}
		Assert.assertTrue(last == 172f);
		Assert.assertTrue(MicroarrayModule.VisualisationMethods.ARRAY_LAYOUT.isApplicableTo(affyMicroarray));

		// SOM
		DataBean somData = manager.createDataBean("som.tsv", new FileInputStream(TestConstants.SOM_CLUSTERED_RESOURCE));
		somData.addTypeTag(MicroarrayModule.TypeTags.SOM_CLUSTERED_EXPRESSION_VALUES);
		somData.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
				
		Assert.assertTrue(MicroarrayModule.VisualisationMethods.SOM.isApplicableTo(somData));
		Table som = somData.queryFeatures("/clusters/som").asTable();
		Assert.assertNotNull(som);
		Assert.assertEquals(som.getColumnCount(), 5); 
		
		// hierarchical clustering
		DataBean hcTree = manager.createDataBean("hs.tre", new FileInputStream(TestConstants.HIERARCHICAL_CLUSTERED_RESOURCE));
		DataBean hcHeatmap = manager.createDataBean("hc.tsv", new FileInputStream(TestConstants.HIERARCHICAL_CLUSTERED_HEATMAP_RESOURCE));
		hcHeatmap.addTypeTag(BasicModule.TypeTags.TABLE_WITH_COLUMN_NAMES);
		hcHeatmap.addTypeTag(MicroarrayModule.TypeTags.NORMALISED_EXPRESSION_VALUES);

		hcTree.addLink(Link.DERIVATION, hcHeatmap);

		Assert.assertTrue(MicroarrayModule.VisualisationMethods.HIERARCHICAL.isApplicableTo(hcTree));
		Table heatmap = hcTree.queryFeatures("/clusters/hierarchical/heatmap").asTable();
		Assert.assertNotNull(heatmap);
		heatmap.nextRow();
		Assert.assertNotNull(heatmap.getValue(" "));
		
		String tree = hcTree.queryFeatures("/clusters/hierarchical/tree").asStrings().iterator().next();
		Assert.assertNotNull(tree);
	}
}
