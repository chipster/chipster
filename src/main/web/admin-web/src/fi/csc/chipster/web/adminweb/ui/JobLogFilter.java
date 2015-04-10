package fi.csc.chipster.web.adminweb.ui;

import java.util.Arrays;
import java.util.Collection;

import com.vaadin.data.hbnutil.filter.ContainerFilter;
import com.vaadin.data.hbnutil.filter.IdContainerFilter;
import com.vaadin.data.hbnutil.filter.StringContainerFilter;
import com.vaadin.event.ShortcutAction;
import com.vaadin.event.ShortcutListener;
import com.vaadin.server.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.NativeSelect;
import com.vaadin.ui.Notification;
import com.vaadin.ui.TextField;

import fi.csc.chipster.web.adminweb.data.DateContainerFilter;
import fi.csc.chipster.web.adminweb.data.JobLogContainer;

public class JobLogFilter extends HorizontalLayout {
	
	
	private static final Collection<String> SEARCH_COLUMNS = Arrays.asList(
			new String[] { 
					JobLogContainer.USERNAME, JobLogContainer.OPERATION, JobLogContainer.COMPHOST, 
					JobLogContainer.START_TIME, JobLogContainer.END_TIME, JobLogContainer.WALLCLOCK_TIME, 
					JobLogContainer.STATUS });
	
	private JobLogView view;
	private TextField searchStringField;
	private NativeSelect columnToSearch;

	private com.vaadin.data.hbnutil.filter.ContainerFilter containerFilter;

	public JobLogFilter(final JobLogView view, String column, String search) {
		this.view = view;

		searchStringField = new TextField();
		if (search != null) {
			searchStringField.setValue(search);
		}
		searchStringField.setDescription("Search for values starting with this string. Question mark (?) is a wildcard for a single character and asterisk (*) for any number of characters.");  
		searchStringField.addShortcutListener(new ShortcutListener("Search", ShortcutAction.KeyCode.ENTER, null) {

			@Override
			public void handleAction(Object sender, Object target) {
				view.update();
			}
		});

		columnToSearch = new NativeSelect();

		Button clearButton = new Button();
		clearButton.setIcon(new ThemeResource("crystal/button_cancel-bw.png"));
		clearButton.setDescription("Remove filter");
		clearButton.addStyleName("search-button");

		for (int i = 0; i < JobLogContainer.NATURAL_COL_ORDER.length; i++) {

			//Do not search from generated columns
			if (SEARCH_COLUMNS.contains(JobLogContainer.NATURAL_COL_ORDER[i])) {
				columnToSearch.addItem(JobLogContainer.NATURAL_COL_ORDER[i]);
				columnToSearch.setItemCaption(JobLogContainer.NATURAL_COL_ORDER[i],
						JobLogContainer.COL_HEADERS_ENGLISH[i]);
			}
		}

		if (column != null) {
			columnToSearch.setValue(column);
		} else {
			columnToSearch.setValue(JobLogContainer.USERNAME);
		}
		columnToSearch.setNullSelectionAllowed(false);

		clearButton.addClickListener(new Button.ClickListener() {
			public void buttonClick(ClickEvent event) {
				getView().clearFilter(JobLogFilter.this);
			}
		});

		addComponent(columnToSearch);
		addComponent(searchStringField);
		addComponent(clearButton);

		addStyleName("search-filter-bg");
		addStyleName("search-filter");
	}

	@Override
	public String toString() {
		return columnToSearch.getValue() + " = *"
				+ (String) searchStringField.getValue() + "*";
	}

	public ContainerFilter getContainerFilter() {
		
		String search = (String) searchStringField.getValue();
						
		if (search == "") {
			containerFilter = null;
			
		} else if (columnToSearch.getValue().equals(JobLogContainer.START_TIME) || columnToSearch.getValue().equals(JobLogContainer.END_TIME)) {

			try {
				containerFilter = new DateContainerFilter(
						columnToSearch.getValue(), search);
			} catch (NumberFormatException e) {
				Notification.show(e.getMessage(), Notification.Type.WARNING_MESSAGE);
				containerFilter = null;
			}

		} else if (columnToSearch.getValue().equals(JobLogContainer.WALLCLOCK_TIME)) {

			try {
			containerFilter = new IdContainerFilter(
					columnToSearch.getValue(), Integer.parseInt(search));
			
			} catch (NumberFormatException e) {
				Notification.show ("Search term "+ columnToSearch.getValue() + " must be numeric", Notification.Type.WARNING_MESSAGE);
				containerFilter = null;
			}
		} else {
			// replace familiar wildcards with SQL LIKE wildcards
			search = search.replace("*", "%");
			search = search.replace("?", "_");
			
			/*
			 * Index is used only when the query is case sensitive and only
			 * prefixes are matched. Users can still override the latter
			 * restriction by adding a wildcard in the beginning, but then the
			 * query will resort to a full table scan.
			 */
			containerFilter = new StringContainerFilter(
					columnToSearch.getValue(), search, false, true);
		}

		return containerFilter;
	}

	public JobLogView getView() {
		return view;
	}

	public void clear() {
		searchStringField.setValue("");
		containerFilter = null;
	}
}