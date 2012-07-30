package fi.csc.microarray.manager.web.ui;

import java.util.Date;

import com.ibm.icu.text.SimpleDateFormat;
import com.vaadin.data.Property;
import com.vaadin.ui.Component;
import com.vaadin.ui.Label;
import com.vaadin.ui.Table;

class DateColumnGenerator implements Table.ColumnGenerator {

	public Component generateCell(Table source, final Object itemId,
			Object columnId) {

		Property prop = source.getItem(itemId).getItemProperty(columnId);
		if (prop != null && prop.getType() != null && prop.getType().equals(Date.class)) {
			
			Date date = (Date) prop.getValue();

			SimpleDateFormat dateFormat = new SimpleDateFormat("yyy-MM-dd HH:mm");
			
			Label label = new Label(dateFormat.format(date));
			return label;
		}
		return null;
	}
}