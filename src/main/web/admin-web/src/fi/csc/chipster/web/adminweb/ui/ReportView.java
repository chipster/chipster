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

	private ChipsterAdminUI app;

	
	public ReportView(ChipsterAdminUI app) {
		
		super(UPDATE_WAIT);
		
		this.app = app;
					
		this.addComponent(getToolbar());
		
		this.addComponent(super.getProggressIndicator());

		tabSheet = new TabSheet();
				
		tabSheet.setSizeFull();
		tabSheet.addSelectedTabChangeListener(new SelectedTabChangeListener() {

			@Override
			public void selectedTabChange(SelectedTabChangeEvent e) {
				updateData();
			}
		});
		
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
	}
	
	public void updateData() {

		if (dataSource == null) {
			dataSource = new ReportDataSource();
		}

		if (tabSheet.getSelectedTab() == filebrokerLayout) {
			updateFileBrokerData();
		}

		if (tabSheet.getSelectedTab() == compLayout) {
			updateCompData();
		}
		
		if (tabSheet.getSelectedTab() == jobmanagerLayout) {
			updateJobmanagerData();
		}
	}
	
	private void updateCompData() {
		
		// comp						
		super.submitUpdateAndWait(new Runnable() {
			
			@Override
			public void run() {				
				dataSource.updateCompReport(ReportView.this);					
			}			
		});
	}

	private void updateFileBrokerData() {

		// filebroker
		super.submitUpdate(new Runnable() {
			
			@Override
			public void run() {				
				dataSource.updateFilebrokerReport(ReportView.this);
			}			
		});
		
	}
	
	private void updateJobmanagerData() {

		super.submitUpdate(new Runnable() {
			
			@Override
			public void run() {				
				dataSource.updateJobmanagerReport(ReportView.this);
			}			
		});
		
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

			toolbarLayout.addComponent(app.getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}

	public void buttonClick(ClickEvent event) {
		if (super.isRefreshButton(event.getSource())) {			
			updateData();			
		}
	}

	public Label getFilebrokerLabel() {
		return filebrokerLabel;
	}
	
	public VerticalLayout getCompLayout() {
		return compLayout;
	}
	
	public VerticalLayout getJobmanagerLayout() {
		return compLayout;
	}
}
