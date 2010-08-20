package fi.csc.microarray.databeans.biobeans;

import java.io.IOException;
import java.util.ArrayList;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.LinkUtils;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.table.EditableTable;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.Modules;

public class PhenodataTest {

	private DataManager manager; 
	
	@BeforeSuite(alwaysRun = true)
	public void init() throws IOException, IllegalConfigurationException, InstantiationException, IllegalAccessException, ClassNotFoundException {
		DirectoryLayout.initialiseSimpleLayout().getConfiguration();			
		this.manager = new DataManager();
		new Modules("fi.csc.microarray.module.chipster.MicroarrayModule").plugFeatures(manager);
	}

	@Test(groups = {"unit"} )
	public void testPhenodataRetrieval() throws MicroarrayException, IOException {
		DataBean normalised = manager.createDataBean("normalised.tsv");
		DataBean filtered = manager.createDataBean("filtered.tsv");
		DataBean filtered2 = manager.createDataBean("filtered2.tsv");
		DataBean phenodata = manager.createDataBean("phenodata.tsv");
		
		filtered.addLink(Link.DERIVATION, normalised);
		filtered2.addLink(Link.DERIVATION, normalised);
		phenodata.addLink(Link.ANNOTATION, normalised);
		
		Assert.assertEquals(phenodata, LinkUtils.retrieveInherited(normalised, Link.ANNOTATION));
		Assert.assertEquals(phenodata, LinkUtils.retrieveInherited(filtered, Link.ANNOTATION));
		Assert.assertEquals(2, LinkUtils.retrieveOutputSet(filtered).length);
	}
	
	
	@Test(groups = {"unit"} )
	public void testPhenodataGeneration() throws MicroarrayException, IOException {

		DataBean normalised1 = manager.createDataBean("normalised.tsv");
		DataBean phenodata1 = manager.createDataBean("phenodata.tsv");
		DataBean filtered = manager.createDataBean("filtered.tsv");
		filtered.addLink(Link.DERIVATION, normalised1);
		phenodata1.addLink(Link.ANNOTATION, normalised1);

		DataBean normalised2 = manager.createDataBean("normalised.tsv");
		DataBean phenodata2 = manager.createDataBean("phenodata.tsv");
		phenodata2.addLink(Link.ANNOTATION, normalised2);
				
		ArrayList<String> samples = new ArrayList<String>();
		ArrayList<String> originals = new ArrayList<String>();
		ArrayList<String> group = new ArrayList<String>();
		ArrayList<String> training = new ArrayList<String>();
		ArrayList<String> chiptypes = new ArrayList<String>();

		samples.add("microarray1.cel");
		originals.add("affy_example1.cel");
		group.add("1");
		training.add("0");
		chiptypes.add("test3");

		samples.add("microarray2.cel");
		originals.add("affy_example2.cel");
		group.add("2");
		training.add("0");
		chiptypes.add("test3");

		EditableTable matrix = new EditableTable(); 
		matrix.addColumn("sample", samples);
		matrix.addColumn("original_name", originals);
		matrix.addColumn("group", group);
		matrix.addColumn("training", training);
		matrix.addColumn("chiptype", chiptypes);

		// FIXME use data manager
		throw new IllegalStateException("Fix these to use DataManager");
//		// write first phenodata out
//		OutputStream out = phenodata1.getContentOutputStreamAndLockDataBean();
//		matrix.writeTo(out);
//		phenodata1.closeContentOutputStreamAndUnlockDataBean(out);
//		
//		// make one group value empty and remove original names
//		group.set(0, "");
//		matrix.removeColumn("group");
//		matrix.addColumn("group", group);
//		matrix.removeColumn("original_name");
//		
//		// write second phenodata out
//		out = phenodata2.getContentOutputStreamAndLockDataBean();
//		matrix.writeTo(out);
//		phenodata2.closeContentOutputStreamAndUnlockDataBean(out);
//		
//		// check phenodata		
//		//   1. complete phenodataset
//		Assert.assertTrue(phenodata1.queryFeatures("/phenodata").exists());
//		Assert.assertTrue(phenodata1.queryFeatures("/phenodata/is-complete").exists());
//		Assert.assertEquals(phenodata1.queryFeatures("/phenodata/describe/microarray1.cel").asString(), "affy_example1.cel");
//		//   2. incomplete phenodataset
//		Assert.assertTrue(phenodata2.queryFeatures("/phenodata").exists());
//		Assert.assertFalse(phenodata2.queryFeatures("/phenodata/is-complete").exists());
//		Assert.assertEquals(phenodata2.queryFeatures("/phenodata/describe/microarray1.cel").asString(), "microarray1.cel");
//		// 3. datasets with links to phenodata
//		Assert.assertTrue(normalised1.queryFeatures("/phenodata/linked").exists());
//		Assert.assertTrue(filtered.queryFeatures("/phenodata/linked/is-complete").exists());
//		Assert.assertEquals(normalised1.queryFeatures("/phenodata/linked/describe/microarray1.cel").asString(), "affy_example1.cel");
//		Assert.assertTrue(normalised2.queryFeatures("/phenodata/linked").exists());

	}
}
