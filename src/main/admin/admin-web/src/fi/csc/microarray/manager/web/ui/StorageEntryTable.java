package fi.csc.microarray.manager.web.ui;

import com.vaadin.ui.Table;

public class StorageEntryTable extends Table {
	public StorageEntryTable(StorageView view) {

		setSelectable(true);
		setImmediate(true);
		addListener(view);
		setNullSelectionAllowed(false);

		setSizeFull();
		
		this.addGeneratedColumn("size", new HumanReadableLongColumnGenerator());
	}
}