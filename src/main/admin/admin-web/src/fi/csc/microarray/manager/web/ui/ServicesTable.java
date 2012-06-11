package fi.csc.microarray.manager.web.ui;

import com.vaadin.ui.Table;

import fi.csc.microarray.manager.web.data.ServiceContainer;

public class ServicesTable extends Table {
	public ServicesTable(ServicesView view) {

		setSelectable(true);
		setImmediate(true);
		addListener(view);
		setNullSelectionAllowed(false);

		setSizeFull();
        
		setContainerDataSource(view.getDataSource());		
	}
}