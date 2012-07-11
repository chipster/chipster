package fi.csc.microarray.manager.web.ui;

import java.util.Map;

import org.hibernate.Session;

import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Button.ClickListener;
import com.vaadin.ui.HorizontalLayout;

import fi.csc.microarray.manager.web.ChipsterAdminApplication;
import fi.csc.microarray.manager.web.data.StatDataSource;
import fi.csc.microarray.manager.web.hbncontainer.JobLogSessionManager;


public class StatView extends HorizontalLayout implements ClickListener {
	
	private ChipsterAdminApplication app;
	private Button testButton;

	public StatView(ChipsterAdminApplication app) {
		this.app = app;
		
		testButton = new Button("Test");
		testButton.addListener(this);
		
		this.addComponent(testButton);
	}

	public void buttonClick(ClickEvent event) {
		if (event.getButton() == testButton) {
			
			Session session = new JobLogSessionManager(app).getSession();
			
			StatDataSource dataSource = new StatDataSource();
			
//			String resultsString = "";
//			
//			for (JobLogEntryRowCount entry : dataSource.getTopUsers(session)) {
//				
//				resultsString += entry.getStartTime() + " " + entry.getRowCount() + "<br>";
//			}
//			
//			app.getMainWindow().showNotification("Results<br>", resultsString);
			
			String resultsString = "";
			
			for (Map map : dataSource.getJobsByMonth(session)) {
				
				for (Object key : map.keySet()) {
				
				resultsString += key + " " + map.get(key) + "<br>";
				}
			}
			
			app.getMainWindow().showNotification("Results<br>", resultsString);
		}
	}
}
