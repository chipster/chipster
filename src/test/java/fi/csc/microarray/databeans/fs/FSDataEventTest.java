package fi.csc.microarray.databeans.fs;

import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.BeforeTest;
import org.testng.annotations.Test;

import fi.csc.microarray.ClientContextUtil;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.databeans.ContentChangedEvent;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataItemCreatedEvent;
import fi.csc.microarray.databeans.DataItemRemovedEvent;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.LinksChangedEvent;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.ModuleManager;

public class FSDataEventTest implements DataChangeListener {
	
	@BeforeTest(groups = {"unit"} )
	public void init() throws Exception {
		ClientContextUtil.setupClientContext();
	}

	@Test(groups = {"unit"} )
	public void testEvents() throws IOException, MicroarrayException, InstantiationException, IllegalAccessException, ClassNotFoundException {
		DataManager manager = new DataManager();
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		moduleManager.plugAll(manager, null);

		moduleManager.plugAll(manager, Session.getSession());
		manager.addDataChangeListener(this);
		manager.setEventsEnabled(true);
		
		DataBean bean1 = manager.createDataBean("My bean.txt");
		ClientContextUtil.setupDatabean(bean1);
		assertNoEvent();
		
		DataBean bean2 = manager.createDataBean("My other bean.txt");
		ClientContextUtil.setupDatabean(bean2);
		assertNoEvent();
		
		bean1.addLink(Link.DERIVATION, bean2);
		assertNoEvent();
		
		manager.connectChild(bean1, manager.getRootFolder());
		assertDataCreatedEvent();
		
		manager.connectChild(bean2, manager.getRootFolder());
		assertDataCreatedEvent();
		
		bean1.addLink(Link.ANNOTATION, bean2);
		assertLinksChangedEvent();
		
		bean1.setName("My nice bean.txt");
		assertContentChangedEvent();
		
		manager.disconnectChild(bean2, manager.getRootFolder());
		assertDataRemovedEvent();
		
	}

	private void assertDataRemovedEvent() {
		Assert.assertFalse(gotTwice);
		Assert.assertNotNull(this.lastEvent);
		Assert.assertTrue(lastEvent instanceof DataItemRemovedEvent, "wrong event type: " + lastEvent.getClass().getSimpleName());
		this.lastEvent = null;
	}

	private void assertContentChangedEvent() {
		Assert.assertFalse(gotTwice);
		Assert.assertNotNull(this.lastEvent);
		Assert.assertTrue(lastEvent instanceof ContentChangedEvent, "wrong event type: " + lastEvent.getClass().getSimpleName());
		this.lastEvent = null;
	}

	private void assertLinksChangedEvent() {
		Assert.assertFalse(gotTwice);
		Assert.assertNotNull(this.lastEvent);
		Assert.assertTrue(lastEvent instanceof LinksChangedEvent, "wrong event type: " + lastEvent.getClass().getSimpleName());
		this.lastEvent = null;
	}

	private void assertDataCreatedEvent() {
		Assert.assertFalse(gotTwice);
		Assert.assertNotNull(this.lastEvent);
		Assert.assertTrue(lastEvent instanceof DataItemCreatedEvent, "wrong event type: " + lastEvent.getClass().getSimpleName());
		this.lastEvent = null;
	}

	private void assertNoEvent() {
		Assert.assertNull(this.lastEvent);
	}

	DataChangeEvent lastEvent = null;
	boolean gotTwice = false;
	
	public void dataChanged(DataChangeEvent event) {
		if (this.lastEvent != null) {
			gotTwice = true; 
		} else {
			this.lastEvent = event;
		}		
	}
}
