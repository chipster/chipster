package fi.csc.microarray.manager.web;

import java.io.IOException;
import java.util.Arrays;

import org.apache.log4j.Logger;

import com.vaadin.Application;
import com.vaadin.ui.Component;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Window;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.manager.Manager;
import fi.csc.microarray.manager.web.ui.JobLogView;
import fi.csc.microarray.manager.web.ui.JobsView;
import fi.csc.microarray.manager.web.ui.ServicesView;
import fi.csc.microarray.manager.web.ui.StatView;
import fi.csc.microarray.manager.web.ui.StorageView;

public class ChipsterAdminApplication extends Application {

	// configuration file path
	//private final String configURL = "http://chipster-devel.csc.fi:8031/chipster-config.xml";
	//private final String configURL = "http://chipster.csc.fi/chipster-config.xml";
	//private final String configURL = "http://chipster.csc.fi/beta/chipster-config.xml";

//	{	
//		try {
//			if (!DirectoryLayout.isInitialised()) {
//				DirectoryLayout.initialiseSimpleLayout(configURL).getConfiguration();
//			}
//		} catch (IOException e) {
//			e.printStackTrace();
//		} catch (IllegalConfigurationException e) {
//			e.printStackTrace();
//		}
//	}

	private final Logger logger = Logger.getLogger(Manager.class);

	private HorizontalLayout horizontalSplit;

	private NavigationMenu navigationLayout;;

	private ServicesView serviceView;
	private StorageView storageView;
	private JobsView jobsView;
	private JobLogView jobLogView;
	private StatView statView;

	private VerticalLayout emptyView;

	private HorizontalLayout toolbarLayout;
	
	public ChipsterAdminApplication() {
		
		// initialise dir and logging
		try {
			//FIXME should obey manager config file given in manager command line parameter
			if (!DirectoryLayout.isInitialised()) {				
				DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"manager"}));
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IllegalConfigurationException e) {
			e.printStackTrace();
		}
		
		//
	}

	private VerticalLayout getServicesView() {
		if (serviceView == null) {

			serviceView = new ServicesView(this);
			try {
				serviceView.loadData();
			} catch (InstantiationException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		return serviceView;
	}
	
	private Component getStorageView() {
		if (storageView == null) {

			storageView = new StorageView(this);
			try {
				storageView.loadData();
			} catch (InstantiationException e) {
				e.printStackTrace();
			} catch (IllegalAccessException e) {
				e.printStackTrace();
			}
		}
		return storageView;
	}

	private JobLogView getJobLogView() {
		if (jobLogView == null) {

			jobLogView = new JobLogView(this);
		}
		return jobLogView;
	}
	
	private JobsView getJobsView() {
		if (jobsView == null) {

			jobsView = new JobsView(this);
		}
		return jobsView;
	}
	
	private StatView getStatView() {
		if (statView == null) {

			statView = new StatView(this);
		}
		return statView;
	}


	@Override
	public void init() {

		
		getJobLogView().init();

		buildMainLayout();
	}



	private void buildMainLayout() {
		setMainWindow(new Window("Chipster admin"));

		horizontalSplit = new HorizontalLayout();
		horizontalSplit.setSizeFull();

		getMainWindow().setContent(horizontalSplit);
		getMainWindow().getContent().setHeight(100, Component.UNITS_PERCENTAGE);
		showEmtpyView();

		setTheme("admin");
	}

	private Component getNavigation() {

		if (navigationLayout == null) {
			navigationLayout = new NavigationMenu(this);
			navigationLayout.showDefaultView();
		}
		return navigationLayout;
	}

	private void setMainComponent(Component c) {

		horizontalSplit.removeAllComponents();

		horizontalSplit.addComponent(getNavigation());
		horizontalSplit.addComponent(c);
		horizontalSplit.setExpandRatio(c, 1);
	}

	protected void showServicesView() {
		setMainComponent(getServicesView());
	}

	public void showJobHistoryView() {
		setMainComponent(getJobLogView());
	}


	public void showJobsView() {
		setMainComponent(getJobsView());
	}


	public void showStorageView() {
		setMainComponent(getStorageView());
	}
	
	public void showStatView() {
		setMainComponent(getStatView());
	}

	public void showEmtpyView() {
		if (emptyView == null) {
			emptyView = new VerticalLayout();
			
			emptyView.addComponent(getToolbar());
			
			emptyView.setSizeFull();
			emptyView.setStyleName("empty-view");
		}
		setMainComponent(emptyView);
	}
	
	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);
			
			toolbarLayout.addComponent(getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setHeight(40, Component.UNITS_PIXELS);
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}
	
	public Component getTitle() {
		Label label = new Label("Chipster admin");
		label.addStyleName("title");
		label.setWidth(250, Component.UNITS_PIXELS);
		return label;
	}
}
