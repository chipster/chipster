package fi.csc.microarray.manager.web.ui;

import com.vaadin.data.Property;
import com.vaadin.terminal.ThemeResource;
import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Component;
import com.vaadin.ui.Table;
import com.vaadin.ui.themes.BaseTheme;

import fi.csc.microarray.manager.web.data.JobLogContainer;

public class JobLogTable extends Table {

	private JobLogView view;

	public JobLogTable(JobLogView view) {

		this.view = view;

		this.setSizeFull();

		this.setWidth("100%");
		this.setSelectable(true);
		this.setImmediate(true);

		this.addGeneratedColumn(JobLogContainer.OUTPUT_LINK, new OutputLinkColumnGenerator());
		this.addGeneratedColumn(JobLogContainer.ERROR_LINK, new ErrorLinkColumnGenerator());
	}

	class OutputLinkColumnGenerator implements Table.ColumnGenerator {

		public Component generateCell(Table source, final Object itemId,
				Object columnId) {

			Button link = new Button("Output");
			link.setStyleName(BaseTheme.BUTTON_LINK);
			link.setDescription("Show job output");

			link.addListener(new Button.ClickListener() {

				public void buttonClick(ClickEvent event) {

					select(itemId);
					view.showOutput(itemId);
				}
			});

			return link;
		}
	}

	class ErrorLinkColumnGenerator implements Table.ColumnGenerator {

		public Component generateCell(Table source, final Object itemId,
				Object columnId) {

			Property prop = source.getItem(itemId).getItemProperty(JobLogContainer.ERROR_MESSAGE);
			if (prop != null && prop.getType() != null && prop.getType().equals(String.class)) {

				String errorMessage = (String) prop.getValue();

				if (errorMessage != null) {

					Button link = new Button();
					link.setIcon(new ThemeResource("../admin/crystal/agt_update_critical.png"));
					link.setStyleName(BaseTheme.BUTTON_LINK);
					link.setDescription("Show job error message");

					link.addListener(new Button.ClickListener() {

						public void buttonClick(ClickEvent event) {

							select(itemId);
							view.showErrorOutput(itemId);
						}
					});

					return link;
				}
			}

			return null;
		}
	}
}
