package fi.csc.microarray.manager.web.ui;

import com.vaadin.ui.Table;

import fi.csc.microarray.manager.web.data.StorageAggregateContainer;

public class StorageAggregateTable extends Table {
	public StorageAggregateTable(StorageView view) {

		setSelectable(true);
		setImmediate(true);
		addListener(view);
		setNullSelectionAllowed(false);
		
		this.setSizeFull();
		
		this.addGeneratedColumn(StorageAggregateContainer.SIZE, new HumanReadableLongColumnGenerator());
	}
}