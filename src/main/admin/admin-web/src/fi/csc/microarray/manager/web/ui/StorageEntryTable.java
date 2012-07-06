package fi.csc.microarray.manager.web.ui;

import com.vaadin.ui.Button;
import com.vaadin.ui.Button.ClickEvent;
import com.vaadin.ui.Component;
import com.vaadin.ui.Table;
import com.vaadin.ui.themes.BaseTheme;

import fi.csc.microarray.manager.web.data.StorageEntryContainer;

public class StorageEntryTable extends Table {
	
	private StorageView view;

	public StorageEntryTable(StorageView view) {
		
		this.view = view;

		setSelectable(true);
		setImmediate(true);
		addListener(view);
		setNullSelectionAllowed(false);

		setSizeFull();
		
		this.addGeneratedColumn(StorageEntryContainer.SIZE, new HumanReadableLongColumnGenerator());
		this.addGeneratedColumn(StorageEntryContainer.DELETE_LINK, new DeleteLinkColumnGenerator());
	}
	
	class DeleteLinkColumnGenerator implements Table.ColumnGenerator {

	    public Component generateCell(Table source, final Object itemId,
	            Object columnId) {
	    	
//	        Property prop = source.getItem(itemId).getItemProperty(columnId);
//	        if (prop != null && prop.getType() != null && prop.getType().equals(Long.class)) {
	        	
	    		Button link = new Button("Delete");
	    		link.setStyleName(BaseTheme.BUTTON_LINK);
	    		
	    		link.addListener(new Button.ClickListener() {

					public void buttonClick(ClickEvent event) {
						
						view.delete(itemId);
					}
	    		});
	    		
	            return link;
//	        }
//	        return null;
	    }
	}
}