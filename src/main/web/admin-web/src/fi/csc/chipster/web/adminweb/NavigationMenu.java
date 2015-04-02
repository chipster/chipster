package fi.csc.chipster.web.adminweb;

import org.apache.log4j.Logger;

import com.vaadin.server.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.Component;
import com.vaadin.ui.Label;
import com.vaadin.ui.VerticalLayout;

public class NavigationMenu extends VerticalLayout implements ClickListener {
	
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(NavigationMenu.class);

	private Button servicesButton = new Button("Services");
	private Button storageButton = new Button("Storage");
	private Button jobsButton = new Button("Jobs");
	private Button jobLogButton = new Button("Job log");
	private Button statButton = new Button("Statistics");
	private Button reportButton = new Button("Maintenance");
	private ChipsterAdminUI app;
	
	private final ThemeResource servicesIcon = new ThemeResource("crystal/ksirtet.png");
	private final ThemeResource servicesIconBw = new ThemeResource("crystal/ksirtet-bw.png");
	private final ThemeResource storageIcon = new ThemeResource("crystal/volume-manager.png");
	private final ThemeResource storageIconBw = new ThemeResource("crystal/volume-manager-bw.png");
	private final ThemeResource jobsIcon = new ThemeResource("crystal/clock.png");
	private final ThemeResource jobsIconBw = new ThemeResource("crystal/clock-bw.png");
	private final ThemeResource jobLogIcon = new ThemeResource("crystal/vcalendar.png");
	private final ThemeResource jobLogIconBw = new ThemeResource("crystal/vcalendar-bw.png");
	private final ThemeResource statIcon = new ThemeResource("crystal/kchart-edited.png");
	private final ThemeResource statIconBw = new ThemeResource("crystal/kchart-edited-bw.png");
	private final ThemeResource reportIcon = new ThemeResource("crystal/service-manager.png");
	private final ThemeResource reportIconBw = new ThemeResource("crystal/service-manager-bw.png");

	public NavigationMenu(ChipsterAdminUI app) {
		this.app = app;

		servicesButton.addClickListener(this);
		addComponent(servicesButton);

		storageButton.addClickListener(this);
		addComponent(storageButton);

		jobsButton.addClickListener(this);
		addComponent(jobsButton);

		jobLogButton.addClickListener(this);
		addComponent(jobLogButton);
		
		statButton.addClickListener(this);
		addComponent(statButton);
		
		reportButton.addClickListener(this);
		addComponent(reportButton);

		setMargin(true);
		setSpacing(true);

		Component placeHolder = new Label(" ");
		addComponent(placeHolder);
		setExpandRatio(placeHolder, 1);

		setHeight("100%");
		setWidth(130, Unit.PIXELS);

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
			reportButton.setIcon(reportIconBw);
			
		} else if (view == storageButton) {
			app.showStorageView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIcon);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIconBw);
			reportButton.setIcon(reportIconBw);
			
		} else if (view == jobsButton) {
			app.showJobsView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIcon);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIconBw);
			reportButton.setIcon(reportIconBw);
			
		} else if (view == jobLogButton) {
			app.showJobLogView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIcon);
			statButton.setIcon(statIconBw);
			reportButton.setIcon(reportIconBw);
			
		} else if (view == statButton) {
			app.showStatView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIcon);
			reportButton.setIcon(reportIconBw);
			
		} else if (view == reportButton) {
			app.showReportView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIconBw);
			reportButton.setIcon(reportIcon);
			
		} else {
			app.showEmtpyView();
			
			servicesButton.setIcon(servicesIconBw);
			storageButton.setIcon(storageIconBw);
			jobsButton.setIcon(jobsIconBw);
			jobLogButton.setIcon(jobLogIconBw);
			statButton.setIcon(statIconBw);
			reportButton.setIcon(reportIconBw);
		}
	}
}
