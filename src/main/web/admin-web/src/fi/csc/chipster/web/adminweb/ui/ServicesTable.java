package fi.csc.chipster.web.adminweb.ui;

import com.vaadin.ui.Table;

import fi.csc.chipster.web.adminweb.data.ServiceContainer;
import fi.csc.chipster.web.adminweb.data.ServiceEntry;

public class ServicesTable extends Table {
	public ServicesTable(ServicesView view) {

		setSelectable(true);
		setNullSelectionAllowed(false);

		setSizeFull();

		setCellStyleGenerator(new CellStyleGenerator() {

			public String getStyle(Table source, Object itemId, Object propertyId) {

				if (itemId instanceof ServiceEntry) {
					ServiceEntry service = (ServiceEntry)itemId;

					if (ServiceContainer.STATUS.equals(propertyId)) {

						if ("UP".equals(service.getStatus().toString())) {
							return "success";
						} else if ("UNKNOWN".equals(service.getStatus().toString())) {
							return "unknown";
						} else {
							return "failure";
						}
					} else if (ServiceContainer.NAME.equals(propertyId)) {

						String name = service.getName().toString();

						if ("authenticator".equals(name) || 
								"comp".equals(name) ||
								"filebroker".equals(name) ||
								"manager".equals(name) ||
								"jobmanager".equals(name) ||
								"client".equals(name)) {
							return name;
						}
					}
				}
				return null;
			}
		});
	}
}