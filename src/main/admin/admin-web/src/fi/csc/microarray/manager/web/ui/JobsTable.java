package fi.csc.microarray.manager.web.ui;

import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Component;
import com.vaadin.ui.Table;
import com.vaadin.ui.themes.BaseTheme;

import fi.csc.microarray.manager.web.data.JobsContainer;

public class JobsTable extends Table {

	private JobsView view;

	public JobsTable(JobsView view) {

		this.view = view;

		this.setSizeFull();

		this.setWidth("100%");
		this.setSelectable(true);
		this.setImmediate(true);

		this.addGeneratedColumn(JobsContainer.CANCEL_LINK, new CancelLinkColumnGenerator());
	}

	class CancelLinkColumnGenerator implements Table.ColumnGenerator {

		public Component generateCell(Table source, final Object itemId,
				Object columnId) {

			Button link = new Button("Cancel");
			link.setStyleName(BaseTheme.BUTTON_LINK);
			link.setDescription("Cancel running job");

			link.addListener(new Button.ClickListener() {

				public void buttonClick(ClickEvent event) {

					select(itemId);
					view.cancel(itemId);
				}
			});

			return link;
		}
	}
}
