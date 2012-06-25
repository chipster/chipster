package fi.csc.microarray.manager.web.ui;

import com.vaadin.data.hbnutil.ContainerFilter;
import com.vaadin.data.hbnutil.StringContainerFilter;
import com.vaadin.event.ShortcutAction;
import com.vaadin.event.ShortcutListener;
import com.vaadin.terminal.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.HorizontalLayout;
import com.vaadin.ui.NativeSelect;
import com.vaadin.ui.TextField;

public class JobLogSearch extends HorizontalLayout {

	private JobLogView view;
	private TextField searchStringField;
	private NativeSelect columnToSearch;

	private StringContainerFilter containerFilter;


	public JobLogSearch(final JobLogView view) {
		this.view = view;

		searchStringField = new TextField();
		searchStringField.addShortcutListener(new ShortcutListener("Search", ShortcutAction.KeyCode.ENTER, null) {

			@Override
			public void handleAction(Object sender, Object target) {
				view.performSearch();
			}
		});


		searchStringField.addStyleName("search-filter-component");
		columnToSearch = new NativeSelect();
		columnToSearch.addStyleName("search-filter-component");

		Button clearButton = new Button();
		clearButton.setIcon(new ThemeResource("crystal/button_cancel-bw.png"));
		clearButton.setDescription("Remove search");
		clearButton.addStyleName("search-button");

		for (int i = 0; i < JobLogView.NATURAL_COL_ORDER.length; i++) {
			columnToSearch.addItem(JobLogView.NATURAL_COL_ORDER[i]);
			columnToSearch.setItemCaption(JobLogView.NATURAL_COL_ORDER[i],
					JobLogView.COL_HEADERS_ENGLISH[i]);
		}

		columnToSearch.setValue("username");
		columnToSearch.setNullSelectionAllowed(false);

		clearButton.addListener(new Button.ClickListener() {
			public void buttonClick(ClickEvent event) {
				getView().clearSearch(JobLogSearch.this);
			}
		});

		addComponent(columnToSearch);
		addComponent(searchStringField);
		addComponent(clearButton);

		setStyleName("search-filter");
	}

	@Override
	public String toString() {
		return columnToSearch.getValue() + " = *"
				+ (String) searchStringField.getValue() + "*";
	}

	public ContainerFilter getContainerFilter() {

		containerFilter = new StringContainerFilter(
				columnToSearch.getValue(), (String) searchStringField.getValue(), true, false);

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