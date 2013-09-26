package fi.csc.chipster.web.adminweb.ui;

import com.vaadin.data.Property;
import com.vaadin.ui.Component;
import com.vaadin.ui.Label;
import com.vaadin.ui.Table;

import fi.csc.chipster.web.adminweb.util.StringUtils;

class HumanReadableLongColumnGenerator implements Table.ColumnGenerator {

    public Component generateCell(Table source, Object itemId,
            Object columnId) {
    	
        Property<?> prop = source.getItem(itemId).getItemProperty(columnId);
        if (prop != null && prop.getType() != null && prop.getType().equals(Long.class)) {
        	
        	Long longSize = (Long)prop.getValue();
        	String stringSize = StringUtils.getHumanReadable(longSize);
        	
            Label label = new Label(stringSize);

            return label;
        }
        return null;
    }
}