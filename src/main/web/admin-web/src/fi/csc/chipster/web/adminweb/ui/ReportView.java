package fi.csc.chipster.web.adminweb.ui;

import java.io.IOException;

import javax.jms.JMSException;

import com.vaadin.shared.ui.label.ContentMode;
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
import fi.csc.chipster.web.adminweb.data.StorageAdminAPI;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;


public class ReportView extends AsynchronousView implements ClickListener {
	
	private VerticalLayout filebrokerLayout = new VerticalLayout();
	private Label filebrokerLabel;
	
	private TabSheet tabSheet;
	private int selectedTab;
	
	private HorizontalLayout toolbarLayout;

	private ReportDataSource dataSource;

	private ChipsterAdminUI app;

	private StorageAdminAPI adminEndpoint;
	
	public ReportView(ChipsterAdminUI app) {
		
		this.app = app;
					
		this.addComponent(getToolbar());
		
		this.addComponent(super.getProggressIndicator());

		tabSheet = new TabSheet();
				
		updateData();
		tabSheet.setSizeFull();
		tabSheet.addSelectedTabChangeListener(new SelectedTabChangeListener() {

			@Override
			public void selectedTabChange(SelectedTabChangeEvent e) {
				
				selectedTab = tabSheet.getTabPosition(tabSheet.getTab(tabSheet.getSelectedTab()));
			}
		});
		
        this.addComponent(tabSheet);        
        this.setExpandRatio(tabSheet, 1);
		this.setSizeFull();
	}
	
	private void updateData() {
		
		try {
			if (dataSource == null) {
				adminEndpoint = new StorageAdminAPI();
				dataSource = new ReportDataSource(adminEndpoint);
			}
			filebrokerLabel = new Label("waiting for status report...", ContentMode.PREFORMATTED);
			filebrokerLabel.addStyleName("report-text");
			
			super.submitUpdate(new Runnable() {

				@Override
				public void run() {				
					dataSource.update(ReportView.this);
				}			
			});

			//selectedTab field is updated when new tabs are added, keep the old value
			int lastSelectedTab = this.selectedTab;

			tabSheet.removeAllComponents();
			filebrokerLayout.removeAllComponents();
			filebrokerLayout.addComponent(filebrokerLabel);
			tabSheet.addTab(filebrokerLayout, "Filebroker");        

			tabSheet.setSelectedTab(lastSelectedTab);
			
		} catch (InstantiationException | IllegalAccessException | IOException | IllegalConfigurationException | MicroarrayException | JMSException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
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
}
