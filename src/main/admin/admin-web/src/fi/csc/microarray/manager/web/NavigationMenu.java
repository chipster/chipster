package fi.csc.microarray.manager.web;

import com.vaadin.terminal.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.Component;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.VerticalLayout;

public class NavigationMenu extends VerticalLayout implements ClickListener {

	private Button servicesButton = new Button("Services");
	private Button storageButton = new Button("Storage");
	private Button jobsButton = new Button("Jobs");
	private Button jobLogButton = new Button("Job log");
	private Button statButton = new Button("Statistics");
	private ChipsterAdminApplication app;
	
	private final ThemeResource servicesIcon = new ThemeResource("crystal/service-manager.png");
	private final ThemeResource servicesIconBw = new ThemeResource("crystal/service-manager-bw.png");
	private final ThemeResource storageIcon = new ThemeResource("crystal/volume-manager.png");
	private final ThemeResource storageIconBw = new ThemeResource("crystal/volume-manager-bw.png");
	private final ThemeResource jobsIcon = new ThemeResource("crystal/clock.png");
	private final ThemeResource jobsIconBw = new ThemeResource("crystal/clock-bw.png");
	private final ThemeResource jobLogIcon = new ThemeResource("crystal/vcalendar.png");
	private final ThemeResource jobLogIconBw = new ThemeResource("crystal/vcalendar-bw.png");
	private final ThemeResource statIcon = new ThemeResource("crystal/kchart-edited.png");
	private final ThemeResource statIconBw = new ThemeResource("crystal/kchart-edited-bw.png");

	public NavigationMenu(ChipsterAdminApplication app) {
		this.app = app;

		servicesButton.addListener(this);
		addComponent(servicesButton);

		storageButton.addListener(this);
		addComponent(storageButton);

		jobsButton.addListener(this);
		addComponent(jobsButton);

		jobLogButton.addListener(this);
		addComponent(jobLogButton);
		
		statButton.addListener(this);
		addComponent(statButton);

		setMargin(true);
		setSpacing(true);
		

		//				navigation.setWidth("100%");
		//
		//				Embedded em = new Embedded("", new ThemeResource("../runo/images/logo.png"));
		//				navigation.addComponent(em);
		//				navigation.setComponentAlignment(em, Alignment.MIDDLE_RIGHT);
		//				navigation.setExpandRatio(em, 1);

		Component placeHolder = new Label(" ");
		addComponent(placeHolder);
		setExpandRatio(placeHolder, 1);

		setHeight("100%");
		setWidth(130, HorizontalLayout.UNITS_PIXELS);

		setStyleName("navigation");
	}
	
	public void showDefaultView() {
		setView(null);
	}
	
	public void buttonClick(ClickEvent event) {
		final Button source = event.getButton();
		
		setView(source);
	}
	
	private void setView(Button view) {
		if (view == servicesButton) {
			app.showServicesView();
			
			servicesButton.setIcon(servicesIcon);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIconBw);
			
		} else if (view == storageButton) {
			app.showStorageView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIcon);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIconBw);
			
		} else if (view == jobsButton) {
			app.showJobsView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIcon);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIconBw);
			
		} else if (view == jobLogButton) {
			app.showJobHistoryView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIcon);
			statButton.setIcon(statIconBw);
			
		} else if (view == statButton) {
			app.showStatView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIcon);
			
		} else {
			app.showEmtpyView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIconBw);
		}
	}
}
