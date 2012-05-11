package fi.csc.microarray.databeans;

import java.io.IOException;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.Traversal;
import fi.csc.microarray.exception.MicroarrayException;

public class LinkTest {
	private DataManager manager; 
	
	@BeforeTest(groups = {"unit"} )
	public void init() throws IOException, IllegalConfigurationException {
		DirectoryLayout.uninitialise();
		DirectoryLayout.initialiseSimpleLayout().getConfiguration();			
		this.manager = new DataManager();
	}
	
	@Test(groups = {"unit"} )
	public void testLinks() throws MicroarrayException {
		DataBean bean1 = manager.createLocalTempDataBean("test1");
		DataBean bean2 = manager.createLocalTempDataBean("test2");
		DataBean bean3 = manager.createLocalTempDataBean("test3");
		
		bean1.addLink(Link.ANNOTATION, bean3);
		bean1.addLink(Link.DERIVATION, bean2);
		bean1.addLink(Link.DERIVATION, bean3);
		bean2.addLink(Link.DERIVATION, bean3);
		
		Assert.assertEquals(bean1.getLinkTargets(Link.DERIVATION).size(), 2);
		Assert.assertEquals(bean1.getLinkSources(Link.DERIVATION).size(), 0);
		Assert.assertEquals(bean2.getLinkTargets(Link.DERIVATION).size(), 1);
		Assert.assertEquals(bean2.getLinkSources(Link.DERIVATION).size(), 1);
		Assert.assertEquals(bean3.getLinkTargets(Link.DERIVATION).size(), 0);
		Assert.assertEquals(bean3.getLinkSources(Link.DERIVATION).size(), 2);
		
		bean1.removeLink(Link.DERIVATION, bean3);
		
		Assert.assertEquals(bean1.getLinkTargets(Link.DERIVATION).size(), 1);
		Assert.assertEquals(bean1.getLinkSources(Link.DERIVATION).size(), 0);
		Assert.assertEquals(bean2.getLinkTargets(Link.DERIVATION).size(), 1);
		Assert.assertEquals(bean2.getLinkSources(Link.DERIVATION).size(), 1);
		Assert.assertEquals(bean3.getLinkTargets(Link.DERIVATION).size(), 0);
		Assert.assertEquals(bean3.getLinkSources(Link.DERIVATION).size(), 1);
		
		Assert.assertEquals(bean1.getLinkTargets(Link.ANNOTATION).size(), 1);
		Assert.assertEquals(bean1.getLinkSources(Link.ANNOTATION).size(), 0);
		Assert.assertEquals(bean2.getLinkTargets(Link.ANNOTATION).size(), 0);
		Assert.assertEquals(bean2.getLinkSources(Link.ANNOTATION).size(), 0);
		Assert.assertEquals(bean3.getLinkTargets(Link.ANNOTATION).size(), 0);
		Assert.assertEquals(bean3.getLinkSources(Link.ANNOTATION).size(), 1);
	}
	
	@Test(groups = {"unit"} )
	public void testTraversal() throws MicroarrayException {
		final DataBean bean1 = manager.createLocalTempDataBean("test1");
		final DataBean bean2 = manager.createLocalTempDataBean("test2");
		final DataBean bean3 = manager.createLocalTempDataBean("test3");
		
		bean1.addLink(Link.DERIVATION, bean2);
		bean1.addLink(Link.DERIVATION, bean3);
		bean2.addLink(Link.DERIVATION, bean3);
		
		Assert.assertTrue(bean1.traverseLinks(Link.derivationalTypes(), Traversal.DIRECT).size() == 3);
		Assert.assertTrue(bean1.traverseLinks(Link.derivationalTypes(), Traversal.REVERSED).size() == 1);
		Assert.assertTrue(bean1.traverseLinks(Link.derivationalTypes(), Traversal.BIDIRECTIONAL).size() == 3);
		Assert.assertTrue(bean2.traverseLinks(Link.derivationalTypes(), Traversal.DIRECT).size() == 2);
		Assert.assertTrue(bean2.traverseLinks(Link.derivationalTypes(), Traversal.REVERSED).contains(bean1));
		Assert.assertTrue(bean2.traverseLinks(Link.derivationalTypes(), Traversal.BIDIRECTIONAL).size() == 3);
		
		List<DataBean> selected = bean2.traverseLinks(Link.derivationalTypes(), Traversal.BIDIRECTIONAL, new DataBeanSelector() {
			public boolean shouldSelect(DataBean bean) {
				return bean == bean2;
			}
			public boolean shouldTraverse(DataBean bean) {
				return true;
			}
		});
		
		Assert.assertTrue(selected.size() == 1);
		Assert.assertTrue(selected.contains(bean2));
	}
}
