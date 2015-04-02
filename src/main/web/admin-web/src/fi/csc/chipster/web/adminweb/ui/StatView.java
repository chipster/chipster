package fi.csc.chipster.web.adminweb.ui;

import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.hibernate.Session;
import org.hibernate.exception.GenericJDBCException;

import com.vaadin.data.Property.ValueChangeEvent;
import com.vaadin.data.Property.ValueChangeListener;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.CheckBox;
import com.vaadin.ui.Component;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.Label;
import com.vaadin.ui.TabSheet;
import com.vaadin.ui.TabSheet.SelectedTabChangeEvent;
import com.vaadin.ui.TabSheet.SelectedTabChangeListener;
import com.vaadin.ui.Table;

import fi.csc.chipster.web.adminweb.ChipsterAdminUI;
import fi.csc.chipster.web.adminweb.data.StatDataSource;
import fi.csc.chipster.web.adminweb.hbncontainer.HibernateUtil;


public class StatView extends AsynchronousView implements ClickListener {
	
	private static final int TIMEOUT = 60;
	
	private Table monthlyStats = new Table();
	private Table yearlyStats = new Table();
	private Table topUsers = new Table();
	private Table toolFails = new Table();
	private Table toolUsage = new Table();
	private Table moduleUsage = new Table();
	
	private TabSheet tabSheet;
	private int selectedTab;
	
	private HorizontalLayout toolbarLayout;
	private CheckBox ignoreTestAccounts;

	private Session session;
	private StatDataSource dataSource;
	
	public StatView(ChipsterAdminUI app) {
		
		super(app, TIMEOUT);
					
		this.addComponent(getToolbar());
		this.addComponent(super.getProggressIndicator());

		tabSheet = new TabSheet();
		tabSheet.setSizeFull();
		
        this.addComponent(tabSheet);        
        this.setExpandRatio(tabSheet, 1);
		this.setSizeFull();
		
        tabSheet.addTab(monthlyStats, "Monthly statistics");        
        tabSheet.addTab(yearlyStats, "Yearly statistics");
        tabSheet.addTab(toolUsage, "Tools usage (1 year)");
        tabSheet.addTab(topUsers, "Top users (1 year)");
        tabSheet.addTab(toolFails, "Tool fails (1 year)");
        tabSheet.addTab(moduleUsage, "Module job counts (beta)");
        
		tabSheet.addSelectedTabChangeListener(new SelectedTabChangeListener() {
			@Override
			public void selectedTabChange(SelectedTabChangeEvent e) {				
				update();
			}
		});
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

	private void updateData(final boolean ignoreTestAccounts, final Component selectedTab) {
		
		if (dataSource == null) {
			dataSource = new StatDataSource();
		}
			
		super.submitUpdate(new Runnable() {				
			@Override
			public void run() {

				if (selectedTab == monthlyStats) {
					List<Map<Object, Object>> stats = dataSource.getMonthlyStats(getHibernateSession(), ignoreTestAccounts);					
					setData(stats, monthlyStats, dataSource.getMonthlyStatsColumnOrder());
				}

				if (tabSheet.getSelectedTab() == yearlyStats) {					
					List<Map<Object, Object>> stats = dataSource.getYearlyStats(getHibernateSession(), ignoreTestAccounts);					
					setData(stats, yearlyStats, dataSource.getYearlyStatsColumnOrder());					
				}

				if (tabSheet.getSelectedTab() == topUsers) {
					List<Map<Object, Object>> stats = dataSource.getTopUsers(getHibernateSession(), ignoreTestAccounts);					
					setData(stats, topUsers, dataSource.getTopUsersColumnOrder());								
				}

				if (tabSheet.getSelectedTab() == toolFails) {
					List<Map<Object, Object>> stats = dataSource.getToolFails(getHibernateSession(), ignoreTestAccounts);					
					setData(stats, toolFails, dataSource.getToolFailsColumnOrder());										
				}

				if (tabSheet.getSelectedTab() == toolUsage) {
					List<Map<Object, Object>> stats = dataSource.getToolUsage(getHibernateSession(), ignoreTestAccounts);					
					setData(stats, toolUsage, dataSource.getToolUsageColumnOrder());								
				}

				if (tabSheet.getSelectedTab() == moduleUsage) {
					List<Map<Object, Object>> stats = dataSource.getModuleUsage(getHibernateSession(), ignoreTestAccounts);					
					setData(stats, moduleUsage, dataSource.getModuleUsageColumnOrder());			
				}
			}			
		});
	}

	protected void setData(final List<Map<Object, Object>> stats, final Table table, final Object[] columnOrder) {

		this.updateUI(new Runnable() {
			public void run() {				
				mapListToTable(stats, table);
				table.setVisibleColumns(columnOrder);
			}
		});
	}

	public HorizontalLayout getToolbar() {

		if (toolbarLayout == null) {
			
			toolbarLayout = new HorizontalLayout();		
			toolbarLayout.addComponent(super.createRefreshButton(this));
								
			ignoreTestAccounts = new CheckBox("Ignore test accounts", true);
			ignoreTestAccounts.addStyleName("toolbar-component");
			toolbarLayout.addComponent(ignoreTestAccounts);
						
			ignoreTestAccounts.addValueChangeListener(new ValueChangeListener() {

				@Override
				public void valueChange(ValueChangeEvent arg0) {
					update();
					tabSheet.setSelectedTab(selectedTab);
				}
			});
			
			Label spaceEater = new Label(" ");
			toolbarLayout.addComponent(spaceEater);
			toolbarLayout.setExpandRatio(spaceEater, 1);

			toolbarLayout.addComponent(getApp().getTitle());	
			
			toolbarLayout.setWidth("100%");
			toolbarLayout.setStyleName("toolbar");
		}

		return toolbarLayout;
	}
	
	private void mapListToTable(List<Map<Object, Object>> list, Table table) {				

		table.setSizeFull();
		table.removeAllItems();
		
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
	}

	public void buttonClick(ClickEvent event) {
		if (super.isRefreshButton(event.getSource())) {			
			update();			
		}
	}

	public void update() {
		updateData(ignoreTestAccounts.getValue(), tabSheet.getSelectedTab());
	}
}
