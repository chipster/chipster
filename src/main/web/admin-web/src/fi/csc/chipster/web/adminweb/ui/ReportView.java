package fi.csc.chipster.web.adminweb.ui;

import org.apache.log4j.Logger;

import com.vaadin.shared.ui.label.ContentMode;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TabSheet;
import com.vaadin.ui.TabSheet.SelectedTabChangeEvent;
import com.vaadin.ui.TabSheet.SelectedTabChangeListener;
import com.vaadin.ui.VerticalLayout;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.ReportDataSource;


public class ReportView extends AsynchronousView implements ClickListener {
	
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(ReportView.class);
	
	protected static final int UPDATE_WAIT = 5; // seconds
	private VerticalLayout filebrokerLayout = new VerticalLayout();
	private VerticalLayout compLayout = new VerticalLayout();
	private VerticalLayout jobmanagerLayout = new VerticalLayout();
	private Label filebrokerLabel;
	private Label compLabel;
	private Label jobmanagerLabel;
	
	private TabSheet tabSheet;
	
	private HorizontalLayout toolbarLayout;

	private ReportDataSource dataSource;
	
	public ReportView(ChipsterAdminUI app) {
		
		super(app, UPDATE_WAIT);
					
		this.addComponent(getToolbar());
		
		this.addComponent(super.getProggressIndicator());

		tabSheet = new TabSheet();
		
		dataSource = new ReportDataSource(app.getEndpoint());
				
		tabSheet.setSizeFull();
		
        this.addComponent(tabSheet);        
        this.setExpandRatio(tabSheet, 1);
		this.setSizeFull();
		
		filebrokerLabel = createReportLabel("waiting for status report...");
		compLabel = createReportLabel("waiting for status report...");
		jobmanagerLabel = createReportLabel("waiting for status report...");
		
		filebrokerLayout.addComponent(filebrokerLabel);
		tabSheet.addTab(filebrokerLayout, "Filebroker");
		compLayout.addComponent(compLabel);
		tabSheet.addTab(compLayout, "Comp");
		jobmanagerLayout.addComponent(jobmanagerLabel);
		tabSheet.addTab(jobmanagerLayout, "Jobmanager");
		
		tabSheet.addSelectedTabChangeListener(new SelectedTabChangeListener() {

			@Override
			public void selectedTabChange(SelectedTabChangeEvent e) {
				update();
			}
		});
	}
	
	public void update() {

		if (tabSheet.getSelectedTab() == filebrokerLayout) {
			super.submitUpdate(new Runnable() {			
				@Override
				public void run() {				
					dataSource.updateFilebrokerReport(ReportView.this);
				}			
			});
		}

		if (tabSheet.getSelectedTab() == compLayout) {
			super.submitUpdateAndWait(new Runnable() {
				@Override
				public void run() {				
					dataSource.updateCompReport(ReportView.this, (int) getTimeout());					
				}			
			});
		}
		
		if (tabSheet.getSelectedTab() == jobmanagerLayout) {
			super.submitUpdate(new Runnable() {			
				@Override
				public void run() {				
					dataSource.updateJobmanagerReport(ReportView.this);
				}			
			});
		}
	}

	public Label createReportLabel(String text) {
		Label label = new Label(text, ContentMode.PREFORMATTED);
		label.addStyleName("report-text");
		
		return label;
	}
	
	public Button createReportButton(String text) {
		Button button = new Button(text);
		button.addStyleName("report-button");
		
		return button;
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			toolbarLayout.addComponent(super.createRefreshButton(this));
														
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);

			toolbarLayout.addComponent(getApp().getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}

	public void buttonClick(ClickEvent event) {
		if (super.isRefreshButton(event.getSource())) {			
			update();			
		}
	}

	public Label getFilebrokerLabel() {
		return filebrokerLabel;
	}
	
	public VerticalLayout getCompLayout() {
		return compLayout;
	}
	
	public VerticalLayout getJobmanagerLayout() {
		return jobmanagerLayout;
	}
}
