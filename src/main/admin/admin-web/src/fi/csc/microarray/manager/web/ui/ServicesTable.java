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
								
				if ("status".equals(propertyId) && itemId instanceof Service) {
					Service service = (Service)itemId;
					
					if ("UP".equals(service.getStatus().toString())) {
						return "success";
					} else {
						return "failure";
					}
				}
				return null;
			}
		});
        
		setContainerDataSource(view.getDataSource());		
	}
}