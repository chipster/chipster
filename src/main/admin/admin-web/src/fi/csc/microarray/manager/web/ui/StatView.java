package fi.csc.microarray.manager.web.ui;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.hibernate.Session;
import org.hibernate.exception.GenericJDBCException;

import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.server.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TabSheet;
import com.vaadin.ui.TabSheet.SelectedTabChangeEvent;
import com.vaadin.ui.TabSheet.SelectedTabChangeListener;
import com.vaadin.ui.Table;
import com.vaadin.ui.VerticalLayout;

import fi.csc.microarray.manager.web.ChipsterAdminUI;
import fi.csc.microarray.manager.web.data.StatDataSource;
import fi.csc.microarray.manager.web.hbncontainer.HibernateUtil;


public class StatView extends VerticalLayout implements ClickListener {
	
	private Table monthlyStats;
	private Table yearlyStats;
	private Table topUsers;
	private Table toolFails;
	private Table toolUsage;
	private Table moduleUsage;
	
	private TabSheet tabSheet;
	private int selectedTab;
	
	private HorizontalLayout toolbarLayout;
	private Button refreshButton;
	private CheckBox ignoreTestAccounts;

	private Session session;
	private StatDataSource dataSource;

	private ChipsterAdminUI app;
	
	public StatView(ChipsterAdminUI app) {
		
		this.app = app;
					
		this.addComponent(getToolbar());

		tabSheet = new TabSheet();
		updateData(ignoreTestAccounts.getValue());
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
	
	private Session getHibernateSession() {
		if (session == null) {
			try {
				session = HibernateUtil.getSessionFactory().openSession();
			} catch (GenericJDBCException e) {
				//FIXME Show exception message and hide or disable all database based content
				e.printStackTrace();			
			}
		}
		return session;
	}

	private void updateData(boolean ignoreTestAccounts) {
		
		Session session = getHibernateSession();
		
		if (dataSource == null) {
			dataSource = new StatDataSource();
		}
		
		monthlyStats = mapListToTable(dataSource.getMonthlyStats(session, ignoreTestAccounts));
		yearlyStats = mapListToTable(dataSource.getYearlyStats(session, ignoreTestAccounts));
		topUsers = mapListToTable(dataSource.getTopUsers(session, ignoreTestAccounts));
		toolFails = mapListToTable(dataSource.getToolFails(session, ignoreTestAccounts));
		toolUsage = mapListToTable(dataSource.getToolUsage(session, ignoreTestAccounts));
		//moduleUsageTable;
		
		monthlyStats.setVisibleColumns(dataSource.getMonthlyStatsColumnOrder());
		yearlyStats.setVisibleColumns(dataSource.getYearlyStatsColumnOrder());
		topUsers.setVisibleColumns(dataSource.getTopUsersColumnOrder());
		toolFails.setVisibleColumns(dataSource.getToolFailsColumnOrder());		
		toolUsage.setVisibleColumns(dataSource.getToolUsageColumnOrder());
		//moduleUsageTable;
		
		//selectedTab field is updated when new tabs are added, keep the old value
		int lastSelectedTab = this.selectedTab;
		
		tabSheet.removeAllComponents();
        tabSheet.addTab(monthlyStats, "Monthly statistics");        
        tabSheet.addTab(yearlyStats, "Yearly statistics");
        tabSheet.addTab(topUsers, "Top users (1 year)");
        tabSheet.addTab(toolFails, "Tool fails (1 year)");
        tabSheet.addTab(toolUsage, "Tools usage (1 year)");
        tabSheet.addTab(moduleUsage, "Module job counts");
                
        tabSheet.setSelectedTab(lastSelectedTab);
	}
	
	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();
			
			refreshButton = new Button("Refresh");
			refreshButton.addClickListener((ClickListener)this);
			refreshButton.setIcon(new ThemeResource("../runo/icons/32/reload.png"));
			refreshButton.setEnabled(true);
			toolbarLayout.addComponent(refreshButton);
								
			ignoreTestAccounts = new CheckBox("Ignore test accounts", true);
			ignoreTestAccounts.addStyleName("toolbar-component");
			toolbarLayout.addComponent(ignoreTestAccounts);
						
			ignoreTestAccounts.addValueChangeListener(new ValueChangeListener() {

				@Override
				public void valueChange(ValueChangeEvent arg0) {
					updateData(ignoreTestAccounts.getValue());
					tabSheet.setSelectedTab(selectedTab);
				}
			});
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);

			toolbarLayout.addComponent(app.getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}
	
	private Table mapListToTable(List<Map<Object, Object>> list) {

		Table table = new Table();
		table.setSizeFull();
		
		if (list.size() > 0) {

			List<Object> keyList = new LinkedList<Object>(list.get(0).keySet());


			//Column headers
			for (Object columnHeader : keyList) {
				table.addContainerProperty(columnHeader, String.class, null);
			}
			
			
			//Content rows
			int i = 0;
			for (Map<Object, Object> map : list) {
			
				//We assume that the order of returned values is identical in all these maps 
				//otherwise values are put to wrong columns
				List<String> stringValues = new LinkedList<String>();
				
				for (Object objValue : map.values()) {
					if (objValue != null) {
						stringValues.add(objValue.toString());
					} else {
						stringValues.add("");
					}					
				}
				
				table.addItem(stringValues.toArray(), i++);
			}
		}
		return table;
	}

	public void buttonClick(ClickEvent event) {
		if (event.getButton() == refreshButton) {			
			updateData(ignoreTestAccounts.getValue());			
		}
	}
}
