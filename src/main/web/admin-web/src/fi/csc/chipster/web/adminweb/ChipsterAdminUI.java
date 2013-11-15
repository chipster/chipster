package fi.csc.chipster.web.adminweb;

import com.vaadin.annotations.Theme;
import com.vaadin.server.ClientConnector.DetachListener;
import com.vaadin.server.VaadinRequest;
import com.vaadin.ui.Component;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.UI;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.adminweb.ui.JobLogView;
import fi.csc.chipster.web.adminweb.ui.JobsView;
import fi.csc.chipster.web.adminweb.ui.ServicesView;
import fi.csc.chipster.web.adminweb.ui.StatView;
import fi.csc.chipster.web.adminweb.ui.StorageView;

@SuppressWarnings("serial")
@Theme("admin")
public class ChipsterAdminUI extends UI implements DetachListener {

	private HorizontalLayout horizontalSplit;

	private NavigationMenu navigationLayout;;

	private ServicesView serviceView;
	private StorageView storageView;
	private JobsView jobsView;
	private JobLogView jobLogView;
	private StatView statView;

	private VerticalLayout emptyView;

	private HorizontalLayout toolbarLayout;
	
	private VerticalLayout getServicesView() {
		if (serviceView == null) {
			serviceView = new ServicesView(this);
		}
		serviceView.update();
		return serviceView;
	}
	
	private Component getStorageView() {
		if (storageView == null) {

			storageView = new StorageView(this);
		}
		storageView.update();
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

	private void buildMainLayout() {

		this.getPage().setTitle("Chipster admin");
		horizontalSplit = new HorizontalLayout();
		horizontalSplit.setSizeFull();

		setContent(horizontalSplit);
		horizontalSplit.setHeight(100, Unit.PERCENTAGE);
		showEmtpyView();	
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

		horizontalSplit.addComponent(getNavigation(), 0);
		horizontalSplit.addComponent(c, 1);
		
		horizontalSplit.setExpandRatio(c, 1);
	}

	protected void showServicesView() {
		setMainComponent(getServicesView());
	}

	public void showJobLogView() {
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
			toolbarLayout.setHeight(40, Unit.PIXELS);
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}
	
	public Component getTitle() {
		Label label = new Label("Chipster admin");
		label.addStyleName("title");
		label.setWidth(250, Unit.PIXELS);
		return label;
	}

	@Override
	protected void init(VaadinRequest request) {
		
		buildMainLayout();
		
		addDetachListener(this);
	}

	@Override
	public void detach(DetachEvent event) {
		if (storageView != null) {
			storageView.clean();
		}
	}
}
