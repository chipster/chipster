package fi.csc.microarray.manager.web;

import java.io.IOException;

import com.vaadin.Application;
import com.vaadin.ui.Component;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.VerticalLayout;
import com.vaadin.ui.Window;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.manager.web.ui.JobLogView;
import fi.csc.microarray.manager.web.ui.JobsView;
import fi.csc.microarray.manager.web.ui.ServicesView;
import fi.csc.microarray.manager.web.ui.StorageView;

public class ChipsterAdminApplication extends Application {

	// configuration file path
	//private final String configURL = "http://chipster-devel.csc.fi:8031/chipster-config.xml";
	private final String configURL = "http://chipster.csc.fi/chipster-config.xml";

	{
		try {
			if (!DirectoryLayout.isInitialised()) {
				DirectoryLayout.initialiseSimpleLayout(configURL).getConfiguration();
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IllegalConfigurationException e) {
			e.printStackTrace();
		}
	}


	private HorizontalLayout horizontalSplit;

	private NavigationMenu navigationLayout;;

	private ServicesView serviceView;
	private StorageView storageView;
	private JobsView jobsView;
	private JobLogView jobLogView;

	private HorizontalLayout emptyView;



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


	@Override
	public void init() {

		getJobLogView().init();

		buildMainLayout();
	}



	private void buildMainLayout() {
		setMainWindow(new Window("Chipster admin web"));

		horizontalSplit = new HorizontalLayout();
		horizontalSplit.setSizeFull();

		getMainWindow().setContent(horizontalSplit);
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

	public void showEmtpyView() {
		if (emptyView == null) {
			emptyView = new HorizontalLayout();
			emptyView.setSizeFull();
			emptyView.setStyleName("empty-view");
		}
		setMainComponent(emptyView);
	}
}
