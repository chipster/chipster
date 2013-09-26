package fi.csc.chipster.web.adminweb.ui;

import com.vaadin.ui.Table;

import fi.csc.chipster.web.adminweb.data.StorageAggregateContainer;

public class StorageAggregateTable extends Table {
	public StorageAggregateTable(StorageView view) {

		setSelectable(true);
		setImmediate(true);
		addValueChangeListener(view);
		setNullSelectionAllowed(false);
		
		this.setSizeFull();
		
		this.addGeneratedColumn(StorageAggregateContainer.SIZE, new HumanReadableLongColumnGenerator());
	}
}