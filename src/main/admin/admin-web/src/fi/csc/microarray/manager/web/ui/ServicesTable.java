package fi.csc.microarray.manager.web.ui;

import com.vaadin.ui.Table;

import fi.csc.microarray.manager.web.data.Service;

public class ServicesTable extends Table {
	public ServicesTable(ServicesView view) {

		setSelectable(true);
		setImmediate(true);
		addListener(view);
		setNullSelectionAllowed(false);

		setSizeFull();

		setCellStyleGenerator(new CellStyleGenerator() {

			public String getStyle(Object itemId, Object propertyId) {

				if (itemId instanceof Service) {
					Service service = (Service)itemId;

					if ("status".equals(propertyId)) {

						if ("UP".equals(service.getStatus().toString())) {
							return "success";
						} else {
							return "failure";
						}
					} else if ("name".equals(propertyId)) {

						String name = service.getName().toString();

						if ("authenticator".equals(name) || 
								"analyser".equals(name) ||
								"filebroker".equals(name) ||
								"manager".equals(name) ||
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