package fi.csc.microarray.manager.web.ui;

import com.vaadin.data.Property;
import com.vaadin.ui.Component;
import com.vaadin.ui.Label;
import com.vaadin.ui.Table;

class HumanReadableLongColumnGenerator implements Table.ColumnGenerator {

    public Component generateCell(Table source, Object itemId,
            Object columnId) {
    	
        Property prop = source.getItem(itemId).getItemProperty(columnId);
        if (prop != null && prop.getType() != null && prop.getType().equals(Long.class)) {
        	
        	Long longSize = (Long)prop.getValue();
        	String stringSize;
        	
        	
        	if (longSize >= 1000000000000l) {
        		stringSize = "" + (longSize / 1000000000000l) + " TB";
        		
        	} else if (longSize >= 1000000000) {
        		stringSize = "" + (longSize / 1000000000) + " GB";
        		
        	} else if (longSize >= 1000000) {
        		stringSize = "" + (longSize / 1000000) + " MB";
        		
        	} else {
        		stringSize = "" + (longSize / 1000) + " kB";
        	}
        	
            Label label = new Label(stringSize);

            return label;
        }
        return null;
    }
}